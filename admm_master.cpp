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
    cout << "Usage: admm_master -X MSNameList.txt -Y outputfile -U uid";
    cout << "-S share data directory ";
    cout << "-t tile size [default:" << Data::TileSize << "]" << endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;
}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "X:Y:P:U:h")) != -1) {
        switch (c) {
            case 'X':
                Data::inputFile = optarg;
                break;
            case 'Y':
                Data::outputFile = optarg;
                break;
            case 'P':
                Data::cmdFile = optarg;
                break;
            case 'U':
                uid = optarg;
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

    if (!Data::inputFile) {
        print_help();
        exit(1);
    }
    Data::MPIData iodata;

    iodata.tilesz = Data::TileSize;
    vector <string> msnames;
    vector <string> inputs;
    vector <struct MsIndex> msIndexs;

    int ntasks = 0;
    ifstream infile(Data::inputFile);

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

    fflush(stdout);

    /**********************************************************/

    iodata.Nms = ntasks;
    /**** get info from slaves ***************************************/

    iodata.freqs = new double[iodata.Nms];
    iodata.freq0 = 0.0;
    iodata.N = iodata.M = iodata.totalt = 0;
    int Mo = 0;
    /* use iodata to store the results, also check for consistency of results */

    for (unsigned int cm = 0; cm < inputs.size(); cm++) {
        FILE *ip = 0;
        if ((ip = fopen(const_cast<char *>(inputs[cm].c_str()),"rb")) == 0) {
            fprintf(stderr, "%s, %s: %d: no input file\n",inputs[cm].c_str(), __FILE__, __LINE__);
            return 1;
        }
        Data::IOData data;
        //fprintf(op, "%s %lf %d %d %d %d %d", Data::TableName, iodata.freq0 * 1e-6, iodata.N, M, Mt, iodata.tilesz, iodata.totalt);
        load_iodata(ip, &data);

        cout << "Slave [" << data.msname << "] N=" << data.N << " M/Mt=" << data.M << "/" << data.Mt << " tilesz="
             << data.tilesz << " totaltime=" << data.totalt << " Freq=" << data.freq0 / 1e6<< "Mhz"<< endl;
        msnames.push_back(string(data.msname));
        struct MsIndex msIndex;
        msIndex.cm = cm;
        ms_short_name(data.msname, msIndex.ms);
        msIndexs.push_back(msIndex);

        if (cm == 0) { /* update data */
            iodata.N = data.N;
            Mo = data.M;
            iodata.M = data.Mt;
            iodata.tilesz = data.tilesz;
            iodata.totalt = data.totalt;
        } else { /* compare against others */
            if ((iodata.N != data.N) || (iodata.M != data.M) || (iodata.tilesz != data.tilesz)) {
                cout << "Slave " << data.msname << " parameters do not match  N=" << data.N << " M=" << data.M
                     << " tilesz=" << data.tilesz << endl;
            }
            if (iodata.totalt < data.totalt) {
                /* use max value as total time */
                iodata.totalt = data.totalt;
            }
        }
        iodata.freqs[cm] = data.freq0;
        iodata.freq0 += data.freq0;

        if(ip) {
            fclose(ip);
        }
        Data::freeData(data);
    }

    iodata.freq0 /= (double) iodata.Nms;

    cout << "Reference frequency (MHz)=" << iodata.freq0 * 1.0e-6 << endl;

    /* regularization factor array, size Mx1 one per each hybrid cluster */
    double *arho, *arhoslave;
    if ((arho = (double *) calloc((size_t) iodata.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((arhoslave = (double *) calloc((size_t) Mo, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }

    /* if text file is given, read it and update rho array */
    if (Data::admm_rho_file) {
        read_arho_fromfile(Data::admm_rho_file, iodata.M, arho, Mo, arhoslave);
    } else {
        /* copy common value */
        /* setup regularization factor array */
        for (int p = 0; p < iodata.M; p++) {
            arho[p] = admm_rho;
        }
        for (int p = 0; p < Mo; p++) {
            arhoslave[p] = admm_rho;
        }
    }
    /* interpolation polynomial */
    double *B, *Bi;
    /* Npoly terms, for each frequency, so Npoly x Nms */
    if ((B = (double *) calloc((size_t) Npoly * iodata.Nms, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(B, 0, sizeof(double) * Npoly * iodata.Nms);
    /* pseudoinverse */
    if ((Bi = (double *) calloc((size_t) Npoly * Npoly, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(Bi, 0, sizeof(double) * Npoly * Npoly);

    setup_polynomials(B, Npoly, iodata.Nms, iodata.freqs, iodata.freq0, (Npoly == 1 ? 1 : PolyType));
    /* determine how many iterations are needed */
    int Ntime = (iodata.totalt + iodata.tilesz - 1) / iodata.tilesz;
    /* override if input limit is given */
    if (Nmaxtime > 0 && Ntime > Nmaxtime) {
        Ntime = Nmaxtime;
    }
    cout << "Master total timeslots=" << Ntime << endl;

    /* send array to slaves */
    /* update rho on each slave */
    iodata.Mo = Mo;

    /*------------------------------------share data file-------------------------------------------------------------*/
    char shareFile[256];
    sprintf(shareFile,"%s/v.mastershare",Data::shareDir);

    FILE *sp = 0;
    if ((sp = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    dump_mpidata(sp, &iodata);
    fwrite(&ntasks, sizeof(int), 1, sp);
    for (unsigned int cm = 0; cm < msIndexs.size(); cm++) {
        fwrite(&(msIndexs[cm]), sizeof(struct MsIndex), 1, sp);
    }
    if (sp) {
        fclose(sp);
    }
    /* dump single mpidata*/
    dump_share_mpidata(Data::shareDir, &iodata);
    write_share_XYZ(Data::shareDir, "master", arho, iodata.M, "arho");
    write_share_XYZ(Data::shareDir, "master", B, Npoly * iodata.Nms, "B");
    write_share_XYZ(Data::shareDir, "master", Bi, Npoly * Npoly, "Bi");

    /* ADMM memory */
    double *Z, *Y, *z;
    /* Z: 2Nx2 x Npoly x M */
    /* keep ordered by M (one direction together) */
    if ((Z = (double *) calloc((size_t) iodata.N * 8 * Npoly * iodata.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(Z, 0, sizeof(double) * iodata.N * 8 * Npoly * iodata.M);

    /* z : 2Nx2 x M x Npoly vector, so each block is 8NM */
    if ((z = (double *) calloc((size_t) iodata.N * 8 * Npoly * iodata.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(z, 0, sizeof(double) * iodata.N * 8 * Npoly * iodata.M);
    /* copy of Y+rho J, M times, for each slave */
    /* keep ordered by M (one direction together) */
    if ((Y = (double *) calloc((size_t) iodata.N * 8 * iodata.M * iodata.Nms, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(Y, 0, sizeof(double) * iodata.N * 8 * iodata.M * iodata.Nms);

    /* need a copy to calculate dual residual */
    double *Zold;
    if ((Zold = (double *) calloc((size_t) iodata.N * 8 * Npoly * iodata.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(Zold, 0, sizeof(double) * iodata.N * 8 * Npoly * iodata.M);

    double *fratio;
    if ((fratio = (double *) calloc((size_t) iodata.Nms, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(fratio, 0, sizeof(double) * iodata.Nms);

    write_share_XYZ(Data::shareDir, "master", Z, iodata.N * 8 * Npoly * iodata.M, "Z");
    write_share_XYZ(Data::shareDir, "master", z, iodata.N * 8 * Npoly * iodata.M, "z");
    write_share_XYZ(Data::shareDir, "master", Y, iodata.N * 8 * iodata.M * iodata.Nms, "Y");
    write_share_XYZ(Data::shareDir, "master", Zold, iodata.N * 8 * Npoly * iodata.M, "Zold");
    write_share_XYZ(Data::shareDir, "master", fratio, iodata.Nms, "fratio");
    /*------------------------------------output data file------------------------------------------------------------*/
    FILE *op = 0;
    if ((op = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    dump_mpidata(op, &iodata);
    fwrite(arhoslave, sizeof(double), iodata.Mo, op);

    if (op) {
        fclose(op);
    }
    /*-----------------------------------------------free ------------------------------------------------------------*/
    delete[] iodata.freqs;
    free(B);
    free(arhoslave);
    free(arho);
    free(Z);
    free(Zold);
    free(z);
    free(fratio);
    cout << "Done." << endl;
    return 0;
}

