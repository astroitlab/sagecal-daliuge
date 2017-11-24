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

void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage:" << endl;
    cout << "write_residual_slave -X inputFile -Y outputfile -P params.txt -U uid" << endl;
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
    /*-----------------------------------------------input------------------------------------------------------------*/
    Data::IOData iodata, old_iodata;
    baseline_t *barr;
    Data::LBeam beam;
    int sources_precessed = get_last_iter(uid);
    int start_iter = 0;
    int M, Mt;

    openblas_set_num_threads(1);
    srand(time(0));

    FILE *ip = 0;
    if ((ip = fopen(Data::inputFile,"rb")) == 0) {
        fprintf(stderr, "%s: %d: no input file\n", __FILE__, __LINE__);
        return 1;
    }
    int flag = 3;
    fread(&flag, sizeof(int), 1, ip); /*--skip flag--*/
    cout << "skip flag:" << flag << endl;
    load_iodata(ip,&old_iodata);
    load_share_iodata(Data::shareDir, old_iodata.msname, &iodata);

    cout << "iodata.N/M/Mt/Nms:" << iodata.N << "/" << iodata.M  << "/" << iodata.Mt << "/" << iodata.Nms
         << ", iodata.freq0:" << iodata.freq0/ 1e6<< "Mhz" << endl;
    Data::freeData(old_iodata);

    if (Data::doBeam) {
        load_share_beam(Data::shareDir, iodata.msname, &beam);
    }

    if ((barr = (baseline_t *) calloc((size_t) iodata.Nbase * iodata.tilesz, sizeof(baseline_t))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    load_share_barr(Data::shareDir, iodata.msname, &iodata, barr);

    if (ip) {
        fclose(ip);
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    clus_source_t *carr;

    read_sky_cluster(Data::SkyModel, Data::Clusters, &carr, &M, iodata.freq0, iodata.ra0, iodata.dec0, Data::format);
    if (M <= 0) {
        fprintf(stderr, "%s: %d: no clusters to solve\n", __FILE__, __LINE__);
        exit(1);
    } else {
        printf("%s:Got %d clusters\n", __FILE__, M);
    }
    Mt = 0;
    int ci = 0, ck = 0,cj = 0;

    for (ci = 0; ci < M; ci++) {
        //printf("cluster %d has %d time chunks\n",carr[ci].id,carr[ci].nchunk);
        Mt += carr[ci].nchunk;
    }
    /* update cluster array with correct pointers to parameters */

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
    /*----------------------------------------------------------------------------------------------------------------*/
    double *xbackup = 0;
    if (iodata.Nchan > 1 || Data::whiten) {
        if ((xbackup = (double *) calloc((size_t) iodata.Nbase * 8 * iodata.tilesz, sizeof(double))) == 0) {
            fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
            exit(1);
        }
        read_share_XYZ(Data::shareDir, iodata.msname, xbackup, iodata.Nbase * 8 * iodata.tilesz, "xbackup");
    }

    complex double *coh;

    if ((coh = (complex double *) calloc((size_t)(iodata.M * iodata.Nbase * iodata.tilesz * 4), sizeof(complex double)))==0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    load_share_coh(Data::shareDir, iodata.msname, &iodata, coh);

    double *p, *pinit;
    /* parameters 8*N*M ==> 8*N*Mt */
    if ((p = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    read_share_XYZ(Data::shareDir, iodata.msname, p, iodata.N * 8 * Mt, "p");

    /* backup of default initial values */
    if ((pinit = (double *) calloc((size_t) iodata.N * 8 * Mt, sizeof(double))) == 0) {
        fprintf(stderr, "%s: %d: no free memory\n", __FILE__, __LINE__);
        exit(1);
    }
    read_share_XYZ(Data::shareDir, iodata.msname, pinit, iodata.N * 8 * Mt, "pinit");
    /*----------------------------------------------------------------------------------------------------------------*/
    int tilex = 0;
    double res_0, res_1, res_00, res_01;
    double mean_nu;
    /* previous residual */
    double res_prev = CLM_DBL_MAX;
    double res_ratio = 15.0; /* how much can the residual increase before resetting solutions, set higher than stand alone mode */
    res_0 = res_1 = res_00 = res_01 = 0.0;

    load_share_res(Data::shareDir, iodata.msname, &start_iter, &res_0, &res_1, &res_00, &res_01, &mean_nu, &tilex);


    /*------------------------------------------processing-------------------------------------------------------------*/

    /* write residuals to output */
    if (!doBeam) {
        calculate_residuals_multifreq(iodata.u, iodata.v, iodata.w, p, iodata.xo, iodata.N, iodata.Nbase,
                                      iodata.tilesz, barr, carr, M, iodata.freqs, iodata.Nchan, iodata.deltaf,
                                      iodata.deltat, iodata.dec0, Data::Nt, Data::ccid, Data::rho, Data::phaseOnly);
    } else {
        calculate_residuals_multifreq_withbeam(iodata.u, iodata.v, iodata.w, p, iodata.xo, iodata.N, iodata.Nbase,
                                               iodata.tilesz, barr, carr, M, iodata.freqs, iodata.Nchan,
                                               iodata.deltaf, iodata.deltat, iodata.dec0,
                                               beam.p_ra0, beam.p_dec0, iodata.freq0, beam.sx, beam.sy,
                                               beam.time_utc, beam.Nelem, beam.xx, beam.yy, beam.zz, Data::Nt,
                                               Data::ccid, Data::rho, Data::phaseOnly);
    }
    tilex += iodata.tilesz;
    /* ---------------------------------------------------print solutions to file ------------------------------------*/
    char ms_full_name[256];
    sprintf(ms_full_name,"%s", iodata.msname);
    char short_name[128];
    ms_short_name(ms_full_name, short_name);
    char solutionFile[256];
    sprintf(solutionFile,"%s/%s.solution",Data::solfile, short_name);

    FILE *sfp = 0;
    /* always create default solution file name MS+'.solutions' */
    if ((sfp = fopen(solutionFile, "a+")) == 0) {
        fprintf(stderr, "%s: %d: no file,solutionFile:%s\n", __FILE__, __LINE__,solutionFile);
        return 1;
    }

    /* write additional info to solution file */
    if (Data::solfile && sources_precessed==0) {
        fprintf(sfp, "# solution file created by SAGECal-Daliuge\n");
        fprintf(sfp, "# freq(MHz) bandwidth(MHz) time_interval(min) stations clusters effective_clusters\n");
        fprintf(sfp, "%lf %lf %lf %d %d %d\n", iodata.freq0 * 1e-6, iodata.deltaf * 1e-6,
                (double) iodata.tilesz * iodata.deltat / 60.0, iodata.N, M, Mt);
    }
    if (Data::solfile) {
        for (cj = 0; cj < iodata.N * 8; cj++) {
            fprintf(sfp, "%d ", cj);
            for (ci = M - 1; ci >= 0; ci--) {
                for (ck = 0; ck < carr[ci].nchunk; ck++) {
                    /* print solution */
                    fprintf(sfp, " %e", p[carr[ci].p[ck] + cj]);
                }
            }
            fprintf(sfp, "\n");
        }
    }

    /* do some quality control */
    /* if residual has increased too much, or all are flagged (0 residual)
      or NaN
      reset solutions to original
      initial values : use residual at 1st ADMM */
    /* do not reset if initial residual is 0, because by def final one will be higher */
    if (res_00 != 0.0 && (res_01 == 0.0 || !isfinite(res_01) || res_01 > res_ratio * res_prev)) {
        cout << "Resetting Solution" << endl;
        /* reset solutions so next iteration has default initial values */
        memcpy(p, pinit, (size_t) iodata.N * 8 * Mt * sizeof(double));
        /* also assume iterations have restarted from scratch */
        start_iter = 1;
        /* also forget min residual (otherwise will try to reset it always) */
        if (res_01 != 0.0 && isfinite(res_01)) {
            res_prev = res_01;
        }
    } else if (res_01 < res_prev) { /* only store the min value */
        res_prev = res_01;
    }

    if (solver_mode == SM_OSLM_OSRLM_RLBFGS || solver_mode == SM_RLM_RLBFGS || solver_mode == SM_RTR_OSRLM_RLBFGS ||
        solver_mode == SM_NSD_RLBFGS) {
        if (Data::verbose) {
            cout << "nu=" << mean_nu << endl;
        }
    }
    if (Data::verbose) {
        cout << iodata.msname << ": Timeslot: " << tilex << " residual: initial=" << res_00 << "/" << res_0 << ",final="
             << res_01 << "/" << res_1  << endl;
    }
    dump_share_barr(Data::shareDir,iodata.msname,&iodata,barr);
    write_share_XYZ(Data::shareDir, iodata.msname, p, iodata.N * 8 * Mt, "p");
    dump_share_res(Data::shareDir, iodata.msname, &start_iter, &res_0, &res_1, &res_00, &res_01, &mean_nu,&tilex);
    dump_share_iodata(Data::shareDir,iodata.msname,&iodata);
    if (Data::doBeam) {
        dump_share_beam(Data::shareDir,iodata.msname,&iodata,&beam);
    }
    dump_share_barr(Data::shareDir,iodata.msname,&iodata,barr);
    /*--------------------------------------------output---------------------------------------------------------------*/
    /* if most data are flagged, only send the original Y we got at the beginning */
    /* for initial ADMM iteration, get back Y with common unitary ambiguity */
    FILE *op = 0;
    if ((op = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    fwrite(&flag, sizeof(int), 1, op);
    if (op) {
        fclose(op);
    }
    /*------------------------------------free ----------------------------------------------------------------------*/

    /* free data memory */
    if (!doBeam) {
        Data::freeData(iodata);
    } else {
        Data::freeData(iodata, beam);
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
    free(barr);
    free(p);
    if (iodata.Nchan > 1 || Data::whiten) {
        free(xbackup);
    }
    free(coh);

    if (sfp) {
        fclose(sfp);
    }
    cout << "Done." << endl;
    return 0;
}

