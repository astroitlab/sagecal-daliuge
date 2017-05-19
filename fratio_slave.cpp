#include "data.h"
#include "cmd.h"
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>

#include <sagecal.h>
#include "utils.h"

using namespace std;
using namespace Data;

char *residualFile;
char *uid;
void print_copyright(void) {
    cout<<"SAGECal-DALIUGE 0.0.1 (C) 2016-2017 CNLAB & ICRAR"<<endl;
}

void print_help(void) {
    cout << "Usage:" << endl;
    cout << "fratio_slave -X inputFile -Y outputfile -Z residualFile -P params.txt -U uid" << endl;
    cout << "Report bugs to <wsl@cnlab.net>" << endl;

}

void ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if (ac < 2) {
        print_help();
        exit(0);
    }
    while ((c = getopt(ac, av, "X:Y:U:Z:P:h")) != -1) {
        switch (c) {
            case 'X':
                Data::inputFile = optarg;
                break;
            case 'Z':
                residualFile = optarg;
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
    Data::IOData old_iodata, iodata;
    Data::LBeam beam;
    Data::MPIData mpiData;

    int sources_precessed = get_last_iter(uid);
    openblas_set_num_threads(1);

    if(sources_precessed!=0) {
        if (get_data_flag(Data::inputFile) != 0) {
            char *tmp = Data::inputFile;
            Data::inputFile = residualFile;
            residualFile = tmp;
        }
    }

    FILE *ip = 0;
    if ((ip = fopen(Data::inputFile,"r")) == 0) {
        fprintf(stderr, "%s: %d: no input file\n", __FILE__, __LINE__);
        return 1;
    }
    load_iodata(ip, &old_iodata);
    load_share_iodata(Data::shareDir, old_iodata.msname, &iodata);
    cout << "iodata.N/M/Mt/Nms:" << iodata.N << "/" << iodata.M  << "/" << iodata.Mt << "/" << iodata.Nms
        << ", iodata.freq0:" << iodata.freq0/ 1e6<< "Mhz" << endl;

    if (Data::doBeam) {
        load_share_beam(Data::shareDir,iodata.msname, &beam);
    }

    load_mpidata(ip,&mpiData);

    if (ip) {
        fclose(ip);
    }
    /*----------------------------------------------------------------------------------------------------------------*/
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

    if(sources_precessed!=0) {
        if (!Data::doBeam) {
            Data::loadData(msitr[0]->table(), old_iodata, &old_iodata.fratio);
        } else {
            Data::loadData(msitr[0]->table(), old_iodata, beam, &old_iodata.fratio);
        }
        Data::writeData(msitr[0]->table(), iodata);
        /* advance to next data chunk */
        for (int cm = 0; cm < iodata.Nms; cm++) {
            for(int it=0;it<sources_precessed;it++)
                (*msitr[cm])++;
        }
    }

    if (!Data::doBeam) {
        Data::loadData(msitr[0]->table(), iodata, &iodata.fratio);
    } else {
        Data::loadData(msitr[0]->table(), iodata, beam, &iodata.fratio);
    }
    /*----------------------------------------------------------------------------------------------------------------*/
    /* downweight factor for regularization, depending on amount of data flagged,
       0.0 means all data are flagged */
    iodata.fratio = 1.0 - iodata.fratio;

    if (Data::verbose) {
        cout << iodata.msname << ": downweight ratio (" << iodata.fratio << ") based on flags.------------------------------------------------------------------" << endl;
    }
    /*----------------------------------------------------------------------------------------------------------------*/
    FILE *op = 0;
    if ((op = fopen(Data::outputFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        return 1;
    }
    dump_share_iodata(Data::shareDir,iodata.msname,&iodata);
    if (Data::doBeam) {
        dump_share_beam(Data::shareDir,iodata.msname,&iodata,&beam);
    }
    dump_iodata(op, &iodata);
    dump_mpidata(op, &mpiData);

    if (op) {
        fclose(op);
    }
    /*------------------------------------free ----------------------------------------------------------------------*/
    for (int cm = 0; cm < iodata.Nms; cm++) {
        delete msitr[cm];
        delete msvector[cm];
    }

    if (!doBeam) {
        Data::freeData(iodata);
    } else {
        Data::freeData(iodata, beam);
    }
    Data::freeData(old_iodata);
    delete[] mpiData.freqs;
    cout << "Done." << endl;
    return 0;
}

