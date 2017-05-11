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

#ifndef LMCUT
#define LMCUT 40
#endif

char *uid;
char *masterFile;

void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage:" << endl;
    cout << "update_y_slave -X inputFile -Y outputfile -Z masterFile -P params.txt -U uid" << endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;

}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "X:Y:Z:U:P:h")) != -1) {
        switch (c) {
            case 'X':
                Data::inputFile = optarg;
                break;
            case 'U':
                uid = optarg;
                break;
            case 'Z':
                masterFile = optarg;
                break;
            case 'P':
                Data::cmdFile = optarg;
                break;
            case 'Y':
                Data::outputFile = optarg;
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

    /*-----------------------------------------------input------------------------------------------------------------*/
    Data::IOData old_iodata, iodata;
    Data::LBeam beam;
    Data::MPIData mpiData;
    double *arho;
    openblas_set_num_threads(1);

    if(get_data_flag(Data::inputFile)!=0) {
        char *tmp = Data::inputFile;
        Data::inputFile = masterFile;
        masterFile = tmp;
        tmp = NULL;
    }

    FILE *ip = 0;
    if ((ip = fopen(Data::inputFile,"rb")) == 0) {
        fprintf(stderr, "%s: %d: no input file\n", __FILE__, __LINE__);
        return 1;
    }
    load_iodata(ip, &old_iodata);
    if (ip) {
        fclose(ip);
    }

    FILE *mp = 0;
    if ((mp = fopen(masterFile,"rb")) == 0) {
        fprintf(stderr, "%s: %d: no master file\n", __FILE__, __LINE__);
        return 1;
    }
    load_mpidata(mp,&mpiData);
    if (mp) {
        fclose(mp);
    }

    load_share_iodata(Data::shareDir, old_iodata.msname, &iodata);
    cout << iodata.msname <<",iodata.N/M/Mt/Nms:" << iodata.N << "/" << iodata.M  << "/" << iodata.Mt << "/" << iodata.Nms
         << ", iodata.freq0:" << iodata.freq0/ 1e6<< "Mhz" << endl;
    if (Data::doBeam) {
        load_share_beam(Data::shareDir, iodata.msname, &beam);
    }
    Data::freeData(old_iodata);

    int M = mpiData.Mo;
    int Mt = mpiData.M;

    if ((arho = (double *) calloc((size_t) mpiData.Mo, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    read_share_XYZ(Data::shareDir, iodata.msname, arho, mpiData.Mo, "arho");

    /*----------------------------------------------------------------------------------------------------------------*/
    clus_source_t *carr;
    read_sky_cluster(Data::SkyModel, Data::Clusters, &carr, &M, iodata.freq0, iodata.ra0, iodata.dec0, Data::format);
    if (M <= 0) {
        fprintf(stderr, "%s: %d: no clusters to solve\n", __FILE__, __LINE__);
        exit(1);
    } else {
        printf("%s:Got %d clusters\n", __FILE__, M);
    }
    /* update cluster array with correct pointers to parameters */
    int ci = 0, ck = 0,cj = 0;
    for (ci = 0; ci < M; ci++) {
        if ((carr[ci].p = (int *) calloc((size_t) carr[ci].nchunk, sizeof(int))) == 0) {
            fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
            exit(1);
        }
        for (ck = 0; ck < carr[ci].nchunk; ck++) {
            carr[ci].p[ck] = cj * 8 * iodata.N;
            cj++;
        }
    }

    double *xbackup = 0;
    if (iodata.Nchan > 1 || Data::whiten) {
        if ((xbackup = (double *) calloc((size_t) iodata.Nbase * 8 * iodata.tilesz, sizeof(double))) == 0) {
            fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
            exit(1);
        }
        read_share_XYZ(Data::shareDir, iodata.msname, xbackup, iodata.Nbase * 8 * iodata.tilesz, "xbackup");
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    complex double *coh;

    if ((coh = (complex double *) calloc((size_t)(iodata.M * iodata.Nbase * iodata.tilesz * 4), sizeof(complex double)))==0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    load_share_coh(Data::shareDir, iodata.msname, &iodata, coh);

    /* ADMM memory */
    double *Z, *Y;
    /* Z: (store B_f Z) 2Nx2 x M */
    if ((Z = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    /* Y, 2Nx2 , M times */
    if ((Y = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    /* primal residual J-BZ */
    double *pres;
    if ((pres = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    double *p;
    /* parameters 8*N*M ==> 8*N*Mt */
    if ((p = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }

    read_share_XYZ(Data::shareDir, iodata.msname, Y, iodata.N * 8 * Mt, "Y");
    read_share_XYZ(Data::shareDir, iodata.msname, Z, iodata.N * 8 * Mt, "Z");
    read_share_XYZ(Data::shareDir, iodata.msname, pres, iodata.N * 8 * Mt, "pres");
    read_share_XYZ(Data::shareDir, iodata.msname, p, iodata.N * 8 * Mt, "p");
    /*----------------------------------------------------------------------------------------------------------------*/
    double res_0, res_1, res_00, res_01;
    double mean_nu;
    res_0 = res_1 = res_00 = res_01 = 0.0;
    int start_iter = 0;

    load_share_res(Data::shareDir, iodata.msname, &start_iter, &res_0, &res_1, &res_00, &res_01, &mean_nu);

    Block<int> sort(1);
    sort[0] = MS::TIME; /* note: only sort over TIME for ms iterator to work */
    vector < MSIter * > msitr;
    vector < MeasurementSet * > msvector;

    MeasurementSet *ms = new MeasurementSet(iodata.msname, Table::Update);
    MSIter *mi = new MSIter(*ms, sort, iodata.deltat * (double) iodata.tilesz);
    msitr.push_back(mi);
    msvector.push_back(ms);

    for (int cm = 0; cm < iodata.Nms; cm++) {
        msitr[cm]->origin();
    }

    if (!Data::doBeam) {
        Data::loadData(msitr[0]->table(), iodata, &iodata.fratio);
    } else {
        Data::loadData(msitr[0]->table(), iodata, beam, &iodata.fratio);
    }

    /*------------------------------------------processing-------------------------------------------------------------*/
    int admm = get_last_iter(uid);
    /* update Y_i <= Y_i + rho (J_i-B_i Z)
        since we already have Y_i + rho J_i, only need -rho (B_i Z) */
    ck = 0;
    for (ci = 0; ci < M; ci++) {
        if (arho[ci] > 0.0) {
            my_daxpy(iodata.N * 8 * carr[ci].nchunk, &Z[ck], -arho[ci], &Y[ck]);
        }
        ck += iodata.N * 8 * carr[ci].nchunk;
    }

    /* calculate primal residual J-BZ */
    my_dcopy(iodata.N * 8 * Mt, p, 1, pres, 1);
    my_daxpy(iodata.N * 8 * Mt, Z, -1.0, pres);

    /* primal residual : per one real parameter */
    /* to remove a load of network traffic and screen output, disable this info */
    if (Data::verbose) {
        cout << iodata.msname << ": ADMM : " << admm << " residual: primal="
             << my_dnrm2(iodata.N * 8 * Mt, pres) / sqrt((double) 8 * iodata.N * Mt) << ", initial=" << res_0
             << ", final=" << res_1 << endl;
    }

    write_share_XYZ(Data::shareDir, iodata.msname, Y, iodata.N * 8 * Mt, "Y");
    write_share_XYZ(Data::shareDir, iodata.msname, Z, iodata.N * 8 * Mt, "Z");
    write_share_XYZ(Data::shareDir, iodata.msname, pres, iodata.N * 8 * Mt, "pres");
    write_share_XYZ(Data::shareDir, iodata.msname, p, iodata.N * 8 * Mt, "p");
    dump_share_res(Data::shareDir, iodata.msname, &start_iter, &res_0, &res_1, &res_00, &res_01, &mean_nu);

    dump_share_iodata(Data::shareDir, iodata.msname, &iodata);
    /*--------------------------------------------output---------------------------------------------------------------*/
    /* if most data are flagged, only send the original Y we got at the beginning */
    /* for initial ADMM iteration, get back Y with common unitary ambiguity */
    FILE *op = 0;
    if ((op = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    int flag = 2;
    fwrite(&flag, sizeof(int), 1, op);
    dump_iodata(op, &iodata);
    fwrite(Y, sizeof(double), iodata.N * 8 * Mt, op);
    fwrite(Z, sizeof(double), iodata.N * 8 * Mt, op);

    if (op) {
        fclose(op);
    }
    /*------------------------------------free ----------------------------------------------------------------------*/
    for (int cm = 0; cm < iodata.Nms; cm++) {
        delete msitr[cm];
        delete msvector[cm];
    }

    exinfo_gaussian *exg;
    exinfo_disk *exd;
    exinfo_ring *exr;
    exinfo_shapelet *exs;
    for (ci = 0; ci < M; ci++) {
        free(carr[ci].ll);
        free(carr[ci].mm);
        free(carr[ci].nn);
        free(carr[ci].sI);
        free(carr[ci].p);
        free(carr[ci].ra);
        free(carr[ci].dec);
        for (cj = 0; cj < carr[ci].N; cj++) {
            /* do a proper typecast before freeing */
            switch (carr[ci].stype[cj]) {
                case STYPE_GAUSSIAN:
                    exg = (exinfo_gaussian *) carr[ci].ex[cj];
                    if (exg) free(exg);
                    break;
                case STYPE_DISK:
                    exd = (exinfo_disk *) carr[ci].ex[cj];
                    if (exd) free(exd);
                    break;
                case STYPE_RING:
                    exr = (exinfo_ring *) carr[ci].ex[cj];
                    if (exr) free(exr);
                    break;
                case STYPE_SHAPELET:
                    exs = (exinfo_shapelet *) carr[ci].ex[cj];
                    if (exs) {
                        if (exs->modes) {
                            free(exs->modes);
                        }
                        free(exs);
                    }
                    break;
                default:
                    break;
            }
        }
        free(carr[ci].ex);
        free(carr[ci].stype);
        free(carr[ci].sI0);
        free(carr[ci].f0);
        free(carr[ci].spec_idx);
        free(carr[ci].spec_idx1);
        free(carr[ci].spec_idx2);
    }
    free(carr);
    free(p);
    if (iodata.Nchan > 1 || Data::whiten) {
        free(xbackup);
    }
    free(Z);
    free(Y);
    free(coh);
    free(arho);
    if (!doBeam) {
        Data::freeData(iodata);
    } else {
        Data::freeData(iodata, beam);
    }
    delete[] mpiData.freqs;
    cout << "Done." << endl;
    return 0;
}

