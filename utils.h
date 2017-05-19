//
// Created by wsl on 4/19/17.
//

#ifndef SAGECAL_DALIUGE_UTILS_H
#define SAGECAL_DALIUGE_UTILS_H

#include "data.h"
#include <sagecal.h>


void dump_iodata(FILE *op, Data::IOData *iodata);

void load_iodata(FILE *op, Data::IOData *iodata);

void dump_share_iodata(char *shareDir, char *msname, Data::IOData *iodata);

void load_share_iodata(char *shareDir, char *msname, Data::IOData *iodata);


void dump_beam(FILE *op, Data::IOData *iodata, Data::LBeam *lBeam);

void load_beam(FILE *op, Data::LBeam *lBeam);

void dump_share_beam(char *shareDir, char *msname, Data::IOData *iodata, Data::LBeam *lBeam);

void load_share_beam(char *shareDir, char *msname, Data::LBeam *lBeam);


void dump_barr(FILE *op, Data::IOData *iodata, baseline_t *barr);

void load_barr(FILE *op,Data::IOData *iodata, baseline_t *barr);

void dump_share_barr(char *shareDir, char *msname, Data::IOData *iodata, baseline_t *barr);

void load_share_barr(char *shareDir, char *msname, Data::IOData *iodata, baseline_t *barr);


void dump_coh(FILE *op, Data::IOData *iodata, complex double *coh);

void load_coh(FILE *op, Data::IOData *iodata, complex double *coh);

void dump_share_coh(char *shareDir, char *msname, Data::IOData *iodata, complex double *coh);

void load_share_coh(char *shareDir, char *msname, Data::IOData *iodata, complex double *coh);

void dump_share_res(char *shareDir, char *msname, int *start_iter, double *res_0, double *res_1, double *res_00, double *res_01, double *mean_nu, int *tilex);

void load_share_res(char *shareDir, char *msname, int *start_iter, double *res_0, double *res_1, double *res_00, double *res_01, double *mean_nu, int *tilex);

void s_dump_data(char *oFile,Data::IOData *iodata, Data::LBeam *lBeam, clus_source_t **carr, baseline_t *barr);

void s_load_data(char *iFile,Data::IOData *iodata, Data::LBeam *lBeam, clus_source_t **carr, baseline_t *barr);

void dump_mpidata(FILE *op, Data::MPIData *iodata);

void load_mpidata(FILE *op, Data::MPIData *iodata);

void dump_share_mpidata(char *shareDir, Data::MPIData *iodata);

void load_share_mpidata(char *shareDir, Data::MPIData *iodata);

int get_last_iter(char *uid);

int get_second_last_iter(char *uid);

void ms_short_name(char *msname, char *shortname);

int get_data_flag(char *fileName);

void read_share_XYZ(char *shareDir, char *msname, double *Y, int len, const char *XYZ);

void write_share_XYZ(char *shareDir, char *msname, double *Y, int len, const char *XYZ);

#endif //SAGECAL_DALIUGE_UTILS_H
