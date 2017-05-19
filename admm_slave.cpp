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

char *arhoFile;
char *uid;

void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage:" << endl;
    cout << "admm_slave -X inputFile -Y outputfile -Z arhoFile -P params.txt -U uid" << endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;

}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "X:Y:P:Z:U:h")) != -1) {
        switch (c) {
            case 'X':
                Data::inputFile = optarg;
                break;
            case 'Z':
                arhoFile = optarg;
                break;
            case 'P':
                Data::cmdFile = optarg;
                break;
            case 'Y':
                Data::outputFile = optarg;
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
                case 's':
                    Data::SkyModel = args[i+1];
                    break;
                case 'c':
                    Data::Clusters = args[i+1];
                    break;
                case 'p':
                    Data::solfile = args[i+1];
                    break;
                case 'g':
                    Data::max_iter = atoi(args[i+1]);
                    break;
                case 'a':
                    Data::DoSim = atoi(args[i+1]);
                    if (Data::DoSim < 0) { Data::DoSim = 1; }
                    break;
                case 'b':
                    Data::doChan = atoi(args[i+1]);
                    if (Data::doChan > 1) { Data::doChan = 1; }
                    break;
                case 'B':
                    Data::doBeam = atoi(args[i+1]);
                    if (Data::doBeam > 1) { Data::doBeam = 1; }
                    break;
                case 'S':
                    Data::shareDir = args[i+1];
                    break;
                case 'F':
                    Data::format = atoi(args[i+1]);
                    if (Data::format > 1) { Data::format = 1; }
                    break;
                case 'e':
                    Data::max_emiter = atoi(args[i+1]);
                    break;
                case 'l':
                    Data::max_lbfgs = atoi(args[i+1]);
                    break;
                case 'm':
                    Data::lbfgs_m = atoi(args[i+1]);
                    break;
                case 'j':
                    Data::solver_mode = atoi(args[i+1]);
                    break;
                case 't':
                    Data::TileSize = atoi(args[i+1]);
                    break;
                case 'I':
                    Data::DataField = args[i+1];
                    break;
                case 'O':
                    Data::OutField = args[i+1];
                    break;
                case 'n':
                    Data::Nt = atoi(args[i+1]);
                    break;
                case 'k':
                    Data::ccid = atoi(args[i+1]);
                    break;
                case 'o':
                    Data::rho = atof(args[i+1]);
                    break;
                case 'L':
                    Data::nulow = atof(args[i+1]);
                    break;
                case 'H':
                    Data::nuhigh = atof(args[i+1]);
                    break;
                case 'R':
                    Data::randomize = atoi(args[i+1]);
                    break;
                case 'W':
                    Data::whiten = atoi(args[i+1]);
                    break;
                case 'J':
                    Data::phaseOnly = atoi(args[i+1]);
                    break;
                case 'x':
                    Data::min_uvcut = atof(args[i+1]);
                    break;
                case 'y':
                    Data::max_uvcut = atof(args[i+1]);
                    break;
                case 'z':
                    Data::ignorefile = args[i+1];
                    break;
                case 'v':
                    Data::verbose = 1;
                    break;
                case 'D':
                    Data::DoDiag = atoi(args[i+1]);
                    if (Data::DoDiag < 0) { Data::DoDiag = 0; }
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

    /*---------------------------------------------2 input------------------------------------------------------------*/
    Data::IOData iodata;
    Data::MPIData mpiData;

    if(get_data_flag(Data::inputFile)!=0) {
        char *tmp = Data::inputFile;
        Data::inputFile = arhoFile;
        arhoFile = tmp;
    }

    FILE *ip = 0;
    if ((ip = fopen(Data::inputFile,"rb")) == 0) {
        fprintf(stderr, "%s: %d: no input file\n", __FILE__, __LINE__);
        return 1;
    }
    load_iodata(ip, &iodata);
    cout << "iodata.N/M/Mt:" << iodata.N << "/" << iodata.M  << "/" << iodata.Mt << ", iodata.freq0:" << iodata.freq0/ 1e6<< "Mhz" << endl;
    if (ip) {
        fclose(ip);
    }

    FILE *op = 0;
    if ((op = fopen(arhoFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no arho file\n", __FILE__, __LINE__);
        return 1;
    }
    load_mpidata(op,&mpiData);
    cout << "mpiData.N:" << mpiData.N << ", mpiData.freq0:" << mpiData.freq0/ 1e6<< "Mhz" << endl;

    double *arho, *arho0;
    if ((arho = (double *) calloc((size_t) mpiData.Mo, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    fread(arho, sizeof(double), mpiData.Mo, op);

    if ((arho0 = (double *) calloc((size_t) mpiData.Mo, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memcpy(arho0, arho, (size_t) mpiData.Mo * sizeof(double));

    if (op) {
        fclose(op);
    }
    /*------------------------------------output----------------------------------------------------------------------*/
    FILE *op2 = 0;
    if ((op2 = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    dump_iodata(op2, &iodata);
    dump_mpidata(op2, &mpiData);
    fwrite(arho, sizeof(double), mpiData.Mo, op);

    if (op2) {
        fclose(op2);
    }
    /*--------------------------------------------  share -------------------------------------------------------*/
    write_share_XYZ(Data::shareDir, iodata.msname, arho, iodata.M, "arho");
    write_share_XYZ(Data::shareDir, iodata.msname, arho0, iodata.M, "arho0");

    double res_0, res_1, res_00, res_01,mean_nu;
    int start_iter = 1, tilex = 0;
    res_0 = res_1 = res_00 = res_01 = mean_nu = 0.0;
    dump_share_res(Data::shareDir, iodata.msname, &start_iter, &res_0, &res_1, &res_00, &res_01, &mean_nu, &tilex);
    double *Z, *Y;
    /* Z: (store B_f Z) 2Nx2 x M */
    if ((Z = (double *) calloc((size_t) iodata.N * 8 * mpiData.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(Z, 0, sizeof(double) * iodata.N * 8 * mpiData.M);
    /* Y, 2Nx2 , M times */
    if ((Y = (double *) calloc((size_t) iodata.N * 8 * mpiData.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(Y, 0, sizeof(double) * iodata.N * 8 * mpiData.M);

    double *pres;
    if ((pres = (double *) calloc((size_t) iodata.N * 8 * mpiData.M, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memset(pres, 0, sizeof(double) * iodata.N * 8 * mpiData.M);



    complex double *coh;
    /* coherencies */
    if ((coh = (complex double *) calloc((size_t)(iodata.M * iodata.Nbase * iodata.tilesz * 4), sizeof(complex double)))==0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    double *xbackup = 0;
    if (iodata.Nchan > 1 || Data::whiten) {
        if ((xbackup = (double *) calloc((size_t) iodata.Nbase * 8 * iodata.tilesz, sizeof(double))) == 0) {
            fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
            exit(1);
        }
        memset(xbackup, 0, sizeof(double) * iodata.Nbase * 8 * iodata.tilesz);
        write_share_XYZ(Data::shareDir, iodata.msname, xbackup, iodata.Nbase * 8 * iodata.tilesz, "xbackup");
    }

    write_share_XYZ(Data::shareDir, iodata.msname, Y, iodata.N * 8 * mpiData.M, "Y");
    write_share_XYZ(Data::shareDir, iodata.msname, Z, iodata.N * 8 * mpiData.M, "Z");

    write_share_XYZ(Data::shareDir, iodata.msname, pres, iodata.N * 8 * mpiData.M, "pres");

    dump_share_coh(Data::shareDir,iodata.msname,&iodata,coh);
    /*--------------------------------------------       free  -------------------------------------------------------*/
    free(Y);
    free(Z);
    free(arho);
    free(pres);
    free(coh);
    free(xbackup);
    Data::freeData(iodata);
    delete[] mpiData.freqs;
    cout << "Done." << endl;
    return 0;
}

