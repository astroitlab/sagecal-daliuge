#include "utils.h"
#include <stdio.h>

void dump_share_iodata(char *shareDir, char *msname, Data::IOData *iodata) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.iodata", shareDir, short_name);

    FILE *op = 0;
    if ((op = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }

    dump_iodata(op,iodata);

    if(op)
        fclose(op);
}
void load_share_iodata(char *shareDir, char *msname, Data::IOData *iodata){
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.iodata", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    load_iodata(op,iodata);

    if(op)
        fclose(op);
}
void dump_iodata(FILE *op, Data::IOData *iodata) {
    int len = 0;
    fwrite(&len, sizeof(int), 1, op);/* this is a flag that indicate share-data file 0*/

    fwrite(&(iodata->M), sizeof(int), 1, op);
    fwrite(&(iodata->Mt), sizeof(int), 1, op);

    len = strlen(iodata->msname);
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->msname, sizeof(char), strlen(iodata->msname), op);

    fwrite(&(iodata->N), sizeof(int), 1, op);
    fwrite(&(iodata->Nbase), sizeof(int), 1, op);
    fwrite(&(iodata->tilesz), sizeof(int), 1, op);
    fwrite(&(iodata->Nchan), sizeof(int), 1, op);
    fwrite(&(iodata->Nms), sizeof(int), 1, op);

    len = iodata->Nms;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->NchanMS, sizeof(int), iodata->Nms, op);

    fwrite(&(iodata->deltat), sizeof(double), 1, op);
    fwrite(&(iodata->totalt), sizeof(int), 1, op);
    fwrite(&(iodata->ra0), sizeof(double), 1, op);
    fwrite(&(iodata->dec0), sizeof(double), 1, op);

    len = iodata->Nbase*iodata->tilesz;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->u, sizeof(double), len, op);
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->v, sizeof(double), len, op);
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->w, sizeof(double), len, op);

    len = 8*len;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->x, sizeof(double), len, op);

    len = len*iodata->Nchan;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->xo, sizeof(double), len, op);

    len = iodata->Nchan;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->freqs, sizeof(double), len, op);

    fwrite(&(iodata->freq0), sizeof(double), 1, op);

    len = iodata->Nbase*iodata->tilesz;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->flag, sizeof(double), len, op);

    fwrite(&(iodata->deltaf), sizeof(double), 1, op);
    fwrite(&(iodata->fratio), sizeof(double), 1, op);
}
void load_iodata(FILE *op, Data::IOData *iodata){
    int len = 0;
    fread(&len, sizeof(int), 1, op);/*skip flag*/

    fread(&(iodata->M), sizeof(int), 1, op);
    fread(&(iodata->Mt), sizeof(int), 1, op);

    fread(&len, sizeof(int), 1, op);
    iodata->msname = (char *)malloc(sizeof(char)*(len+1));
    fread(iodata->msname, sizeof(char), len, op);
    iodata->msname[len] = '\0';

    fread(&(iodata->N), sizeof(int), 1, op);
    fread(&(iodata->Nbase), sizeof(int), 1, op);
    fread(&(iodata->tilesz), sizeof(int), 1, op);
    fread(&(iodata->Nchan), sizeof(int), 1, op);
    fread(&(iodata->Nms), sizeof(int), 1, op);

    fread(&len, sizeof(int), 1, op);
    iodata->NchanMS = new int[len];
    fread(iodata->NchanMS, sizeof(int), len, op);

    fread(&(iodata->deltat), sizeof(double), 1, op);
    fread(&(iodata->totalt), sizeof(int), 1, op);
    fread(&(iodata->ra0), sizeof(double), 1, op);
    fread(&(iodata->dec0), sizeof(double), 1, op);

    fread(&len, sizeof(int), 1, op);
    iodata->u = new double[len];
    fread(iodata->u, sizeof(double), len, op);

    fread(&len, sizeof(int), 1, op);
    iodata->v = new double[len];
    fread(iodata->v, sizeof(double), len, op);

    fread(&len, sizeof(int), 1, op);
    iodata->w = new double[len];
    fread(iodata->w, sizeof(double), len, op);

    fread(&len, sizeof(int), 1, op);
    iodata->x = new double[len];
    fread(iodata->x, sizeof(double), len, op);

    fread(&len, sizeof(int), 1, op);
    iodata->xo = new double[len];
    fread(iodata->xo, sizeof(double), len, op);

    fread(&len, sizeof(int), 1, op);
    iodata->freqs = new double[len];
    fread(iodata->freqs, sizeof(double), len, op);

    fread(&(iodata->freq0), sizeof(double), 1, op);

    fread(&len, sizeof(int), 1, op);
    iodata->flag = new double[len];
    fread(iodata->flag, sizeof(double), len, op);

    fread(&(iodata->deltaf), sizeof(double), 1, op);
    fread(&(iodata->fratio), sizeof(double), 1, op);
}

void dump_share_beam(char *shareDir, char *msname, Data::IOData *iodata, Data::LBeam *lBeam){
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.beam", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    dump_beam(op, iodata, lBeam);

    if(op)
        fclose(op);
}
void load_share_beam(char *shareDir, char *msname, Data::LBeam *lBeam){
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.beam", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    load_beam(op,lBeam);
    if(op)
        fclose(op);
}
void dump_beam(FILE *op, Data::IOData *iodata, Data::LBeam *lBeam){

    int len = iodata->tilesz;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(lBeam->time_utc, sizeof(double), len, op);

    len = iodata->N;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(lBeam->Nelem, sizeof(int), len, op);

    fwrite(lBeam->sx, sizeof(double), len, op);

    fwrite(lBeam->sy, sizeof(double), len, op);

    fwrite(lBeam->sz, sizeof(double), len, op);

    fwrite(lBeam->len_xyz, sizeof(int), len, op);

    for(int i=0;i< iodata->N;i++) {
        fwrite(lBeam->xx[i], sizeof(double), lBeam->len_xyz[i], op);
        fwrite(lBeam->yy[i], sizeof(double), lBeam->len_xyz[i], op);
        fwrite(lBeam->zz[i], sizeof(double), lBeam->len_xyz[i], op);
    }

    fwrite(&(lBeam->p_ra0), sizeof(double), 1, op);
    fwrite(&(lBeam->p_dec0), sizeof(double), 1, op);

}
void load_beam(FILE *op, Data::LBeam *lBeam){
    int len = 0;
    fread(&len, sizeof(int), 1, op);
    lBeam->time_utc = new double[len];
    fread(lBeam->time_utc, sizeof(double), len, op);

    fread(&len, sizeof(int), 1, op);
    lBeam->Nelem = new int[len];
    fread(lBeam->Nelem, sizeof(int), len, op);

    lBeam->sx = new double[len];
    fread(lBeam->sx, sizeof(double), len, op);

    lBeam->sy = new double[len];
    fread(lBeam->sy, sizeof(double), len, op);

    lBeam->sz = new double[len];
    fread(lBeam->sz, sizeof(double), len, op);

    lBeam->len_xyz = new int[len];
    fread(lBeam->len_xyz, sizeof(int), len, op);

    lBeam->xx = new double *[len];
    lBeam->yy = new double *[len];
    lBeam->zz = new double *[len];

    for(int i=0; i < len;i++) {
        lBeam->xx[i] = new double[lBeam->len_xyz[i]];
        lBeam->yy[i] = new double[lBeam->len_xyz[i]];
        lBeam->zz[i] = new double[lBeam->len_xyz[i]];
        fread(lBeam->xx[i], sizeof(double), lBeam->len_xyz[i], op);
        fread(lBeam->yy[i], sizeof(double), lBeam->len_xyz[i], op);
        fread(lBeam->zz[i], sizeof(double), lBeam->len_xyz[i], op);
    }

    fread(&(lBeam->p_ra0), sizeof(double), 1, op);
    fread(&(lBeam->p_dec0), sizeof(double), 1, op);
}

void dump_share_barr(char *shareDir, char *msname, Data::IOData *iodata, baseline_t *barr){
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.barr", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    dump_barr(op, iodata, barr);
    if(op)
        fclose(op);
}

void load_share_barr(char *shareDir, char *msname, Data::IOData *iodata, baseline_t *barr){
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.barr", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    load_barr(op, iodata, barr);
    if(op)
        fclose(op);
}
void dump_barr(FILE *op, Data::IOData *iodata, baseline_t *barr){
    for(int i=0;i<iodata->Nbase * iodata->tilesz;i++) {
        fwrite(&(barr[i].sta1), sizeof(int), 1, op);
        fwrite(&(barr[i].sta2), sizeof(int), 1, op);
        fwrite(&(barr[i].flag), sizeof(unsigned char), 1, op);
    }
}

void load_barr(FILE *op, Data::IOData *iodata, baseline_t *barr){
    for(int i=0;i<iodata->Nbase * iodata->tilesz;i++) {
        fread(&(barr[i].sta1), sizeof(int), 1, op);
        fread(&(barr[i].sta2), sizeof(int), 1, op);
        fread(&(barr[i].flag), sizeof(unsigned char), 1, op);
    }
}
void dump_coh(FILE *op, Data::IOData *iodata, complex double *coh) {
    fwrite(coh, sizeof(complex double), iodata->M * iodata->Nbase * iodata->tilesz * 4, op);
}

void load_coh(FILE *op, Data::IOData *iodata, complex double *coh) {
    fread(coh, sizeof(complex double), iodata->M * iodata->Nbase * iodata->tilesz * 4, op);
}

void dump_share_coh(char *shareDir, char *msname, Data::IOData *iodata, complex double *coh) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.coh", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    dump_coh(op, iodata, coh);
    if(op)
        fclose(op);
}

void load_share_coh(char *shareDir, char *msname, Data::IOData *iodata, complex double *coh) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.coh", shareDir, short_name);

    FILE *op = 0;

    if ((op = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    load_coh(op, iodata, coh);
    if(op)
        fclose(op);
}
void dump_share_mpidata(char *shareDir, Data::MPIData *iodata) {
    char shareFile[256];
    sprintf(shareFile,"%s/v.share.mpidata", shareDir);

    FILE *op = 0;

    if ((op = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    dump_mpidata(op, iodata);
    if(op)
        fclose(op);
}
void s_dump_data(char *oFile,Data::IOData *iodata, Data::LBeam *lBeam, clus_source_t **carr, baseline_t *barr) {
    FILE *op = 0;
    if ((op = fopen(oFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        exit(1);
    }
    dump_iodata(op, iodata);
    if(Data::doBeam) {
        dump_beam(op, iodata, lBeam);
    }
    dump_barr(op, iodata, barr);

    if (op) {
        fclose(op);
    }
}

void s_load_data(char *iFile,Data::IOData *iodata, Data::LBeam *lBeam, clus_source_t **carr, baseline_t *barr) {
    FILE *op = 0;
    if ((op = fopen(iFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no output file\n", __FILE__, __LINE__);
        exit(1);
    }
    load_iodata(op, iodata);
    if(Data::doBeam) {
        load_beam(op, lBeam);
    }

    load_barr(op, iodata, barr);
    if (op) {
        fclose(op);
    }
}

void load_share_mpidata(char *shareDir, Data::MPIData *iodata) {
    char shareFile[256];
    sprintf(shareFile,"%s/v.share.mpidata", shareDir);

    FILE *op = 0;

    if ((op = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    load_mpidata(op, iodata);
    if(op)
        fclose(op);
}

void dump_mpidata(FILE *op, Data::MPIData *iodata) {
    int len = 1;
    fwrite(&len, sizeof(int), 1, op);/* this is a flag that indicate master share-data file 1*/

    fwrite(&(iodata->N), sizeof(int), 1, op);
    fwrite(&(iodata->M), sizeof(int), 1, op);
    fwrite(&(iodata->tilesz), sizeof(int), 1, op);
    fwrite(&(iodata->Nms), sizeof(int), 1, op);
    fwrite(&(iodata->totalt), sizeof(int), 1, op);

    len = iodata->Nms;
    fwrite(&len, sizeof(int), 1, op);
    fwrite(iodata->freqs, sizeof(double), len, op);

    fwrite(&(iodata->freq0), sizeof(double), 1, op);
    fwrite(&(iodata->Mo), sizeof(int), 1, op);
}

void load_mpidata(FILE *op, Data::MPIData *iodata) {
    int len = 0;
    fread(&len, sizeof(int), 1, op);/*skip flag*/

    fread(&(iodata->N), sizeof(int), 1, op);
    fread(&(iodata->M), sizeof(int), 1, op);
    fread(&(iodata->tilesz), sizeof(int), 1, op);
    fread(&(iodata->Nms), sizeof(int), 1, op);
    fread(&(iodata->totalt), sizeof(int), 1, op);

    fread(&len, sizeof(int), 1, op);
    iodata->freqs = new double[len];
    fread(iodata->freqs, sizeof(double), len, op);
    fread(&(iodata->freq0), sizeof(double), 1, op);
    fread(&(iodata->Mo), sizeof(int), 1, op);

}

int get_last_iter(char *uid){
    char v_uid[256];
    sprintf(v_uid,"%s", uid);
    char * pch = strtok (v_uid, "/");
    char *idx;
    while(pch!= NULL) {
        idx = pch;
        pch = strtok (NULL, "/");
    }
    return atoi(idx);
}

int get_second_last_iter(char *uid){
    char v_uid[256];
    sprintf(v_uid,"%s", uid);
    char * pch = strtok (v_uid, "/");
    char *idx;
    do {
        pch = strtok (NULL, "/");
        if ( pch != NULL) {
            idx = pch;
            pch = strtok (NULL, "/");
            if ( pch == NULL) {
                break;
            }
        }
    } while ( pch!= NULL);
    return atoi(idx);
}

void ms_short_name(char *msname, char *shortname) {
    char * pch = strtok (msname, "/");
    char *v_shortname = pch;
    while ( pch!= NULL) {
        v_shortname = pch;
        pch = strtok (NULL, "/");
    }
    snprintf(shortname, strlen(v_shortname)+1, "%s", v_shortname);
}

int get_data_flag(char *fileName) {
    FILE *op = 0;
    if ((op = fopen(fileName, "rb")) == 0) {
        fprintf(stderr, "%s: %d: not found file\n", __FILE__, __LINE__);
        return -1;
    }
    int len;
    fread(&len, sizeof(int), 1, op);/* this is a flag that indicate master share-data file 1*/
    if (op) {
        fclose(op);
    }
    return len;
}

void read_share_XYZ(char *shareDir, char *msname, double *Y, int len, const char *XYZ) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.%s", shareDir, short_name, XYZ);

    FILE *sp = 0;

    if ((sp = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        if ((sp = fopen(shareFile, "wb")) == 0) {
            exit(1);
        }
        fclose(sp);
        return;
    }
    fread(Y, sizeof(double),len,sp);
    if (sp) {
        fclose(sp);
    }
}

void write_share_XYZ(char *shareDir, char *msname, double *Y, int len, const char *XYZ) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];

    sprintf(shareFile,"%s/v.%s.share.%s", shareDir, short_name, XYZ);
    FILE *sp = 0;

    if ((sp = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    fwrite(Y, sizeof(double),len,sp);
    if (sp) {
        fclose(sp);
    }
}

void dump_share_res(char *shareDir, char *msname, int *start_iter, double *res_0, double *res_1, double *res_00, double *res_01, double *mean_nu, int *tilex) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.res", shareDir, short_name);

    FILE *sp = 0;

    if ((sp = fopen(shareFile, "wb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    fwrite(start_iter, sizeof(int),1,sp);
    fwrite(res_0, sizeof(double),1,sp);
    fwrite(res_1, sizeof(double),1,sp);
    fwrite(res_00, sizeof(double),1,sp);
    fwrite(res_01, sizeof(double),1,sp);
    fwrite(mean_nu, sizeof(double),1,sp);
    fwrite(tilex, sizeof(int),1,sp);
    if (sp) {
        fclose(sp);
    }
}

void load_share_res(char *shareDir, char *msname, int *start_iter, double *res_0, double *res_1, double *res_00, double *res_01, double *mean_nu, int *tilex) {
    char short_name[64], ms_full_name[256];
    sprintf(ms_full_name,"%s", msname);
    ms_short_name(ms_full_name,short_name);
    char shareFile[300];
    sprintf(shareFile,"%s/v.%s.share.res", shareDir, short_name);

    FILE *sp = 0;

    if ((sp = fopen(shareFile, "rb")) == 0) {
        fprintf(stderr, "%s: %d: no share data file %s\n", __FILE__, __LINE__, shareFile);
        exit(1);
    }
    fread(start_iter, sizeof(int),1,sp);
    fread(res_0, sizeof(double),1,sp);
    fread(res_1, sizeof(double),1,sp);
    fread(res_00, sizeof(double),1,sp);
    fread(res_01, sizeof(double),1,sp);
    fread(mean_nu, sizeof(double),1,sp);
    fread(tilex, sizeof(int),1,sp);
    if (sp) {
        fclose(sp);
    }
}