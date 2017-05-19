#include "data.h"
#include "cmd.h"
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "sagecal.h"
#include "utils.h"

using namespace std;
using namespace Data;
char *uid;

void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage:" << endl;
    cout << "admm_slave -X inputFile -Y outputfile -P params.txt -U uid" << endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;
}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "X:Y:U:P:h")) != -1) {
        switch (c) {
            case 'X':
                Data::inputFile = optarg;
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
    ifstream infile(Data::inputFile);
    /* check if the file exists and readable */
    if (!infile.good()) {
        cout << "Input File " << Data::inputFile << " does not exist." << endl;
        exit(1);
    }

    string buffer;
    if (infile.is_open()) {
        if (infile.good()) {
            std::getline(infile, buffer);
            if (buffer.length() > 0) {
                char table_name[256];
                snprintf(table_name, buffer.length(), "%s", buffer.c_str());
                Data::TableName = table_name;
            }
        }
    }
    if (Data::SkyModel && Clusters && Data::TableName) {
        cout << "Using SkyModel: " << Data::SkyModel << ", Clusters: " << Clusters << ", Data::TableName: " << Data::TableName << endl;
    } else {
        cout << "ERROR,  SkyModel: " << Data::SkyModel << ", Clusters: " << Clusters << ", Data::TableName: " << Data::TableName << endl;
    }
    if (!Data::SkyModel || !Data::Clusters || !Data::TableName) {
        exit(1);
    }
    cout << "Selecting baselines > " << Data::min_uvcut << " and < " << Data::max_uvcut << " wavelengths." << endl;
    if (!DoSim) {
        cout << "Using ";
        if (Data::solver_mode == SM_LM_LBFGS || Data::solver_mode == SM_OSLM_LBFGS || Data::solver_mode == SM_RTR_OSLM_LBFGS ||
                Data::solver_mode == SM_NSD_RLBFGS) {
            cout << "Gaussian noise model for solver." << endl;
        } else {
            cout << "Robust noise model for solver with degrees of freedom [" << Data::nulow << "," << Data::nuhigh << "]." << endl;
        }
    } else {
        cout << "Only doing simulation (with possible correction for cluster id " << Data::ccid << ")." << endl;
    }
}

int main(int argc, char **argv) {
    ParseCmdLine(argc, argv);

    Data::IOData iodata;
    Data::LBeam beam;
    iodata.tilesz=Data::TileSize;
    iodata.deltat=1.0;
    iodata.msname = (char *)malloc(sizeof(char) * (strlen(Data::TableName))+1);
    snprintf(iodata.msname, strlen(Data::TableName)+1, "%s", Data::TableName);
    openblas_set_num_threads(1);

    if (!doBeam) {
        readAuxData(iodata.msname, &iodata);
    } else {
        readAuxData(iodata.msname, &iodata, &beam);
    }

    /**********************************************************/
    int M, Mt, ci, cj, ck;
    double *p, *pinit;
    double **pm;

    clus_source_t *carr;
    baseline_t *barr;
    read_sky_cluster(Data::SkyModel, Data::Clusters, &carr, &M, iodata.freq0, iodata.ra0, iodata.dec0, Data::format);
    /* exit if there are 0 clusters (incorrect sky model/ cluster file)*/
    if (M <= 0) {
        fprintf(stderr, "%s: %d: no clusters to solve\n", __FILE__, __LINE__);
        exit(1);
    } else {
        printf("Got %d clusters\n", M);
    }
    /* array to store baseline->sta1,sta2 map */
    if ((barr = (baseline_t *) calloc((size_t) iodata.Nbase * iodata.tilesz, sizeof(baseline_t))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    generate_baselines(iodata.Nbase, iodata.tilesz, iodata.N, barr, Data::Nt);

    /* calculate actual no of parameters needed,
     this could be > M */
    Mt = 0;
    for (ci = 0; ci < M; ci++) {
        //printf("cluster %d has %d time chunks\n",carr[ci].id,carr[ci].nchunk);
        Mt += carr[ci].nchunk;
    }
    printf("Total effective clusters: %d\n", Mt);

    /* parameters 8*N*M ==> 8*N*Mt */
    if ((p = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    /* update cluster array with correct pointers to parameters */
    cj = 0;
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
    /* pointers to parameters */
    if ((pm = (double **) calloc((size_t) Mt, sizeof(double *))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    /* setup the pointers */
    for (ci = 0; ci < Mt; ci++) {
        pm[ci] = &(p[ci * 8 * iodata.N]);
    }
    /* initilize parameters to [1,0,0,0,0,0,1,0] */
    for (ci = 0; ci < Mt; ci++) {
        for (cj = 0; cj < iodata.N; cj++) {
            pm[ci][8 * cj] = 1.0;
            pm[ci][8 * cj + 6] = 1.0;
        }
    }
    /* backup of default initial values */
    if ((pinit = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    memcpy(pinit, p, (size_t) iodata.N * 8 * Mt * sizeof(double));

    /**********************************************************/
    /* timeinterval in seconds */
    cout << "For " << iodata.tilesz << " samples, solution time interval (s): "
         << iodata.deltat * (double) iodata.tilesz << endl;
    cout << "Freq: " << iodata.freq0 / 1e6 << " MHz, Chan: " << iodata.Nchan << " Bandwidth: " << iodata.deltaf / 1e6
         << " MHz" << ", Nms: " << iodata.Nms << endl;

    iodata.M = M;
    iodata.Mt = Mt;
    /*------------------------------------share data file-------------------------------------------------------------*/
    dump_share_iodata(Data::shareDir, iodata.msname, &iodata);
    if (doBeam) {
        dump_share_beam(Data::shareDir, iodata.msname, &iodata, &beam);
    }
    dump_share_barr(Data::shareDir, iodata.msname, &iodata, barr);
    write_share_XYZ(Data::shareDir, iodata.msname, p, iodata.N * 8 * Mt, "p");
    write_share_XYZ(Data::shareDir, iodata.msname, pinit, iodata.N * 8 * Mt, "pinit");
    /*------------------------------------output----------------------------------------------------------------------*/
    FILE *op = 0;
    if ((op = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    dump_iodata(op, &iodata);
    if (doBeam) {
        dump_beam(op, &iodata, &beam);
    }
    if (op) {
        fclose(op);
    }
    /*------------------------------------free ----------------------------------------------------------------------*/

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
    free(barr);
    free(p);
    free(pinit);
    free(pm);
    /* free data memory */
    if (!doBeam) {
        Data::freeData(iodata);
    } else {
        Data::freeData(iodata, beam);
    }

    cout << "Done." << endl;
    return 0;
}

