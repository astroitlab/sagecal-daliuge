#include "data.h"
#include "cmd.h"
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>

#include "sagecal.h"
#include "utils.h"

using namespace std;
using namespace Data;

char *uid;

void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage: write_z_master -X inputFile -Y outputFile -P params_master.txt -U uid ";
    cout << " -p solutions.txt";
    cout << " -S share directory"<< endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;
}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "U:X:Y:P:h")) != -1) {
        switch (c) {
            case 'U':
                uid = optarg;
                break;
            case 'X':
                Data::inputFile = optarg;
                break;
            case 'Y':
                Data::outputFile = optarg;
                break;
            case 'P':
                Data::cmdFile = optarg;
                break;
            case 'h':
                print_help();
                exit(1);
            default:
                print_help();
                exit(1);
        }
    }

    ifstream infile(cmdFile);
    /* check if the file exists and readable */
    if (!infile.good()) {
        cout << "File " << Data::cmdFile << " does not exist." << endl;
        exit(1);
    }
    int argc=0;
    char **args = new char *[108];
    string buffer;
    if (infile.is_open()) {
        while (infile.good()) {
            std::getline(infile, buffer);
            if (buffer.length() > 0) {
                char *cmdline = new char[buffer.size() + 1];
                strcpy(cmdline, buffer.c_str());
                cmdline[buffer.size()] = '\0';
                char * pch = strtok (cmdline, " ");
                args[argc++] = pch;
                while ( pch!= NULL) {
                    pch = strtok (NULL, " ");
                    args[argc++] = pch;
                }
            }
        }

        for(int i=0;i<argc-2;i=i+2) {
            if( args[i][0]!= '-') {
                i = i - 1;
                continue;
            }
            char c = args[i][1];
            switch (c) {
                case 't':
                    Data::TileSize = atoi(args[i+1]);
                    break;
                case 'p':
                    Data::solfile = args[i+1];
                    break;
                case 'r':
                    admm_rho= atof(args[i+1]);
                    break;
                case 'G':
                    admm_rho_file= args[i+1];
                case 'N':
                    Data::Npoly = atoi(args[i+1]);
                    break;
                case 'S':
                    Data::shareDir = args[i+1];
                    break;
                case 'Q':
                    PolyType= atoi(args[i+1]);
                    break;
                case 'v':
                    Data::verbose = 1;
                    break;
            }
        }
        CheckParams(argc-1, args);
    }

}

void CheckParams(int ac, char **av) {

}

int main(int argc, char **argv) {
    ParseCmdLine(argc, argv);
    int first_iter = get_last_iter(uid);
    Data::MPIData iodata;
    load_share_mpidata(Data::shareDir, &iodata);

    FILE *sfp = 0;
    if (solfile) {
        char solutionFile[256];
        sprintf(solutionFile,"%s/master.solution",Data::solfile);
        if ((sfp = fopen(solutionFile, "a+")) == 0) {
            fprintf(stderr, "%s: %d: no file\n", __FILE__, __LINE__);
            return 1;
        }
    }
    /* write additional info to solution file */
    if (solfile && first_iter == 0) {
        fprintf(sfp, "# solution file (Z) created by SAGECal-Daliuge\n");
        fprintf(sfp, "# reference_freq(MHz) polynomial_order stations clusters effective_clusters\n");
        fprintf(sfp, "%lf %d %d %d %d\n", iodata.freq0 * 1e-6, Npoly, iodata.N, iodata.Mo, iodata.M);
        fclose(sfp);
        
        FILE *op = 0;
        if ((op = fopen(Data::outputFile, "wb")) == 0) {
            fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        } else {
            int flag = 2;
            fwrite(&flag, sizeof(int), 1, op);
            fclose(op);
        }

        delete[] iodata.freqs;
        cout << "Done." << endl;
        return 0;
    }

    /*------------------------------------share ----------------------------------------------------------------------*/
    double *Z;
    /* Z: 2Nx2 x Npoly x M */
    /* keep ordered by M (one direction together) */
    if ((Z = (double *) calloc((size_t) iodata.N * 8 * Npoly * iodata.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }

    read_share_XYZ(Data::shareDir, (char *)"master", Z, iodata.N * 8 * Npoly * iodata.M, "Z");
    /*----------------------------------------processing--------------------------------------------------------------*/

    /* write Z to solution file, same format as J, but we have Npoly times more
       values per timeslot per column */
    if (solfile) {
        for (int p = 0; p < iodata.N * 8 * Npoly; p++) {
            fprintf(sfp, "%d ", p);
            for (int pp = 0; pp < iodata.M; pp++) {
                fprintf(sfp, " %e", Z[pp * iodata.N * 8 * Npoly + p]);
            }
            fprintf(sfp, "\n");
        }
        fclose(sfp);
    }
    int resetcount = 0, msgcode = 0;
    for (int cm = 0; cm < iodata.Nms; cm++) {
        /* to-do : how to fetch msgcode from slave*/
        if (msgcode == CTRL_RESET) {
            resetcount++;
        }
    }
    if (resetcount > iodata.Nms / 2) {
        /* if most slaves have reset, print a warning only */
        //memset(Z,0,sizeof(double)*(size_t)iodata.N*8*Npoly*iodata.M);
        cout << "Resetting Global Solution" << endl;
    }

    /*-------------------------------------------output---------------- ----------------------------------------------*/
    FILE *op = 0;
    if ((op = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
    } else {
        int flag = 2;
        fwrite(&flag, sizeof(int), 1, op);
        fclose(op);
    }
    /*-------------------------------------------free  ---------------------------------------------------------------*/
    delete[] iodata.freqs;
    free(Z);
    cout << "Done." << endl;
    return 0;
}

