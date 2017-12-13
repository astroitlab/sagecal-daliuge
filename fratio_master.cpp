//
// Created by wsl on 4/19/17.
//
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

char *slaveFile = NULL;
char *solutionFile = NULL;
char *uid = NULL;

void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage:" << endl;
    cout << "fratio_master -X inputFile -Y outputfile -P params.txt -U uid" << endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;
}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "X:Y:Z:V:U:P:h")) != -1) {
        switch (c) {
            case 'X':
                Data::inputFile = optarg;
                break;
            case 'Y':
                Data::outputFile = optarg;
                break;
            case 'U':
                uid = optarg;
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

    Data::MPIData mpiData;
    int ntasks = 0;
    vector <struct MsIndex> msIndexs;

    /*------------------------------------share data file-------------------------------------------------------------*/
    char shareFile[256];
    sprintf(shareFile,"%s/v.mastershare",Data::shareDir);

    FILE *sp = 0;
    if ((sp = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    load_mpidata(sp, &mpiData);
    fread(&ntasks, sizeof(int), 1, sp);

    for (unsigned int cm = 0; cm <ntasks; cm++) {
        struct MsIndex msIndex;
        fread(&msIndex, sizeof(struct MsIndex), 1, sp);
        msIndexs.push_back(msIndex);
    }
    for (unsigned int cm = 0; cm < msIndexs.size(); cm++) {
        cout << "msindex:," <<msIndexs[cm].cm << "," << msIndexs[cm].ms << endl;
    }
    if (sp) {
        fclose(sp);
    }
    /*------------------------------------input ----------------------------------------------------------------------*/
    double *fratio;
    if ((fratio = (double *) calloc((size_t) mpiData.Nms, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    double *B, *Bi;
    /* Npoly terms, for each frequency, so Npoly x Nms */
    if ((B = (double *) calloc((size_t) Npoly * mpiData.Nms, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    /* pseudoinverse */
    if ((Bi = (double *) calloc((size_t) Npoly * Npoly, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    read_share_XYZ(Data::shareDir, "master", fratio, mpiData.Nms, "fratio");
    read_share_XYZ(Data::shareDir, "master", B, Npoly * mpiData.Nms, "B");
    read_share_XYZ(Data::shareDir, "master", Bi, Npoly * Npoly, "Bi");

    ifstream infile(Data::inputFile);
    vector <string> inputs;
    /* check if the file exists and readable */
    if(!infile.good()) {
        cout <<"File "<<Data::inputFile<<" does not exist."<<endl;
        exit(1);
    }
    string buffer;
    if (infile.is_open()) {
        while(infile.good()) {
            std::getline(infile, buffer);
            if (buffer.length()>0) {
                inputs.push_back(buffer);
                ntasks = ntasks + 1;
            }
        }
    }
    for (unsigned int cm = 0; cm < inputs.size(); cm++) {
        FILE *ip = 0;
        if ((ip = fopen(const_cast<char *>(inputs[cm].c_str()), "rb")) == 0) {
            fprintf(stderr, "%s, %s: %d: no input file\n", inputs[cm].c_str(), __FILE__, __LINE__);
            return 1;
        }
        Data::IOData data;
        //fprintf(op, "%s %lf %d %d %d %d %d", Data::TableName, iodata.freq0 * 1e-6, iodata.N, M, Mt, iodata.tilesz, iodata.totalt);
        load_iodata(ip, &data);
        char short_name[128], ms_full_name[256];
        sprintf(ms_full_name,"%s", data.msname);
        ms_short_name(ms_full_name, short_name);

        for (unsigned int sm = 0; sm < msIndexs.size(); sm++) {
            if(strcmp(short_name, msIndexs[sm].ms)==0) {
                cout << "find:" << msIndexs[sm].ms << ", idx:" << sm << endl;
                fratio[msIndexs[sm].cm] = data.fratio;
                break;
            }
        }
        Data::freeData(data);
        if(ip) {
            fclose(ip);
        }
    }
    /*---------------------------------------find_prod_inverse--------------------------------------------------------*/
    /* interpolation polynomial */

    find_prod_inverse(B, Bi, Npoly, mpiData.Nms, fratio);

    /*----------------------------------------share output------------------------------------------------------------*/
    write_share_XYZ(Data::shareDir, "master", fratio, mpiData.Nms, "fratio");
    write_share_XYZ(Data::shareDir, "master", B, Npoly * mpiData.Nms, "B");
    write_share_XYZ(Data::shareDir, "master", Bi, Npoly * Npoly, "Bi");

    /*-------------------------------------------output---------------------------------------------------------------*/
    FILE *op = 0;
    if ((op = fopen(Data::outputFile,"wb")) == 0) {
        fprintf(stderr, "%s: %d: no input file\n", __FILE__, __LINE__);
        return 1;
    }
    dump_mpidata(op, &mpiData);
    if (op) {
        fclose(op);
    }
    /*---------------------------------------------   free  ----------------------------------------------------------*/
    free(B);
    free(Bi);
    free(fratio);
    delete[] mpiData.freqs;

    cout << "Done." << endl;
    return 0;
}

