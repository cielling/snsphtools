/*
   PURPOSE:
    to convert the units in an SDF file
    This version should also work for large files, since it reads in the
    sdf file line-by-line.

   COMPILE:
    with Makefile: 
    make ARCH=<arch-value> [CC=<c-compiler>] PROGS=getyield [-f Makefile]
    This routine needs some libraries from the tree code and SDF routines,
    so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
    once for serial use (i.e. without PAROS flag).

   RUN:
    convertunits <in-file.sdf> <convert.ctl> <out-file.sdf>

   NOTES:
    the <convert.ctl> should contain the conversion factors for
     mass
     length
     time
    in that order. Factors should be given as the ratio old-unit/new-unit.

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <SDF.h>


typedef enum SDF_type_enum SDF_type_t;

typedef union {
    int i;
    float f;
    double d;
} datum_t;

float masscf, lengthcf, timecf;

static void writeinit(FILE *fp);
static void writescalars(SDF *sdfp, FILE *fp);
static void writestructs(SDF *sdfp, FILE *fp);
static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp, FILE **fpctl);


int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;
    FILE *fpctl = NULL;

    initargs(argc, argv, &sdfp, &fp, &fpctl);

    writeinit(fp);

    fscanf(fpctl, "%e %e %e", &masscf, &lengthcf, &timecf);
    printf("read conversion factors:\nmass: %e\nlength: %e\ntime: %e\n",
           masscf, lengthcf, timecf);
    fclose(fpctl);

    writescalars(sdfp, fp);
    writestructs(sdfp, fp);

    fclose(fp);
    SDFclose(sdfp);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp, FILE **fpctl)
{
    char input;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s SDFfile ctl-file outfile \n", argv[0]);
        exit(1);
    }

    *sdfp = SDFopen(NULL, argv[1]);
    if (*sdfp == NULL) {
        fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
        exit(2);
    }

    if (access(argv[3], F_OK) == 0) {
        fprintf(stderr, "%s: %s exists; overwrite (y/n)? ", argv[0],
            argv[2]);
        input = getc(stdin);
        if ((input != 'y') && (input != 'Y'))
            exit(3);
    }
    *fp = fopen(argv[3], "w");
    if (*fp == NULL) {
        fprintf(stderr, "%s: %s\n", argv[3], strerror(errno));
        exit(errno);
    }

    *fpctl = fopen(argv[2], "r");
    if (*fp == NULL) {
        fprintf(stderr, "%s: %s\n", argv[2], strerror(errno));
        exit(errno);
    }

}

static void writeinit(FILE *fp)
{
    fprintf(fp, "# SDF\n");
    fprintf(fp, "parameter byteorder = %#x;\n", SDFcpubyteorder());
}

static void writescalars(SDF *sdfp, FILE *fp)
{
    int i, nvecs, nrecs;
    int flag = 0;
    char **vecs;
    SDF_type_t type;
    datum_t datum;

    nvecs = SDFnvecs(sdfp);/*figure out number of lines in the header, basically-CIE*/
    vecs = SDFvecnames(sdfp);/*get the names of the variables/parameters in the header-CIE*/

    /*go through all of them individually-CIE*/
    for (i = 0; i < nvecs; ++i) {
        /*figure out if scalar(=1) or array(!=1)?-CIE*/
/*
        nrecs = SDFnrecs(vecs[i], sdfp);
        if (nrecs != 1) continue;
*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
        if (flag) continue;

        type = SDFtype(vecs[i], sdfp);
        /* Doesn't handle string */
        /* Also, SDFgetfloat was happy to read a double scalar and
           convert it for me; that probably works backwards too.  I
           don't think there's an equivalent for the SDFrdvecs family
           though, so "read; convert type; write" is the general
           path. */
        /*read in header file, line by line, with the appropriate function-CIE*/
        switch (type) {
        case SDF_INT:
        /*SDFget*:sdfp=sdf file; vecs[i]= variable name; datum.*=holds value of that variable-CIE*/
            SDFgetint(sdfp, vecs[i], &(datum.i));
            break;
        case SDF_FLOAT:
            SDFgetfloat(sdfp, vecs[i], &(datum.f));
            break;
        case SDF_DOUBLE:
            SDFgetdouble(sdfp, vecs[i], &(datum.d));
            break;
        default:
            fprintf(stderr, "%s: type not supported\n", vecs[i]);
/*             exit(-1); */
        }

        if( strncmp(vecs[i], "Gnewt", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f*lengthcf*lengthcf*lengthcf/(masscf*timecf*timecf);

        } else if ( strncmp(vecs[i], "tpos", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "tvel", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "t_wind", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * timecf;

        } else if ( strncmp(vecs[i], "R0", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * lengthcf;
        
        } else if ( strncmp(vecs[i], "bndry_x", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_y", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_z", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_r", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * lengthcf;
        
        } else if ( strncmp(vecs[i], "bndry_vx", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_vy", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_vz", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * lengthcf/ timecf;
        
        } else if ( strncmp(vecs[i], "bndry_px", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_py", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_pz", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * lengthcf*masscf/timecf;
        
        } else if ( strncmp(vecs[i], "bndry_lx", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_ly", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "bndry_lz", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * lengthcf*masscf/timecf;
        
        } else if ( strncmp(vecs[i], "bndry_mass", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f * masscf;
        
        } else if ( strncmp(vecs[i], "ke", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "pe", strlen(vecs[i])) == 0 ||
                    strncmp(vecs[i], "te", strlen(vecs[i])) == 0 ) {
            datum.d = datum.d * masscf*lengthcf*lengthcf/(timecf*timecf);

        } else if( strncmp(vecs[i], "massCF", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f /masscf;

        } else if( strncmp(vecs[i], "lengthCF", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f /lengthcf;

        } else if( strncmp(vecs[i], "timeCF", strlen(vecs[i])) == 0 ) {
            datum.f = datum.f /timecf;

        }
    
        
        /*write header file, line by line, as the appropriate data type-CIE*/
        switch (type) {
        case SDF_INT:
            fprintf(fp, "int %s = %d;\n", vecs[i], datum.i);
            break;
        case SDF_FLOAT:
            fprintf(fp, "float %s = %.7g;\n", vecs[i], datum.f);
            break;
        case SDF_DOUBLE:
            fprintf(fp, "double %s = %.16g;\n", vecs[i], datum.d);
            break;
        default:
            fprintf(stderr, "%s: type not supported\n", vecs[i]);
/*             exit(-1); */
        }
    }
}

static void writestructs(SDF *sdfp, FILE *fp)
{
    FILE *frp = NULL;
    int i, j, nvecs, nmembers;
    int len, counter, flag = 0;
    int nrecs;
    int *strides, *nobjs, *starts, *inoffsets;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0;
    void *btab;
    void **addrs;
    float tmp;
    double tmp_d;

    frp = fopen("log.out", "w");
    if(!frp) printf("error opening log file!\n");

    SDFgetint(sdfp, "npart", &nrecs);
    fprintf(frp,"%d particles\n", nrecs);

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
        if (flag) ++nmembers;
    }
    fprintf(frp,"nmembers = %d\n",nmembers);


    /*malloc memory space for the respective features of the struct-CIE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    strides = (int *)malloc(nmembers * sizeof(int));
    nobjs = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));

    /*one by one, go through all the fields in the column, i.e. members of the struct?-CIE*/
    flag = 0;
    for (i = 0, nmembers = 0, stride = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
        if (flag) {
            members[nmembers] = vecs[i];
            nobjs[nmembers] = 1;/*nobjs[0] is the number of particles. are all elements
                                      in this nobjs array the same then??-CIE:yes, but can be
                                      different*/
            starts[nmembers] = 0;  /* Not correct in parallel; use
                                      NobjInitial */
            types[nmembers] = SDFtype(members[nmembers], sdfp);
            inoffsets[nmembers] = stride;
            stride += SDFtype_sizes[types[nmembers]];

            ++nmembers;
        }
    }

    btab = (void *)malloc(stride);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
    for (i = 0; i < nmembers; ++i) {
        addrs[i] = (char *)btab + inoffsets[i];
        strides[i] = stride;
    }


    /*print the struct declaration part from the header-CIE*/
    fprintf(fp, "struct {\n");
    for (i = 0; i < nmembers; ++i) {
        switch (types[i]) {
        case SDF_INT:
            fprintf(fp, "\tint %s;\n", members[i]);
            break;
        case SDF_FLOAT:
            fprintf(fp, "\tfloat %s;\n", members[i]);
            break;
        case SDF_DOUBLE:
            fprintf(fp, "\tdouble %s;\n", members[i]);
            break;
        default:
            fprintf(stderr, "%s: type not supported\n", members[i]);
            exit(-1);
        }
    }
    fprintf(fp, "}[%d];\n", nrecs);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    for (j = 0; j < nrecs; ++j) {

        for( i = 0; i < nmembers; i++) starts[i] = j;
        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
                 strides);

        /* do conversion */
        for (i = 0; i < nmembers; i++) {

            if( strncmp(members[i],"x",strlen(members[i])) == 0 ||
                strncmp(members[i],"y",strlen(members[i])) == 0 ||
                strncmp(members[i],"z",strlen(members[i])) == 0 ) {
                
                tmp_d = *(double *)(btab + inoffsets[i]);
                tmp_d = tmp_d * lengthcf;
                memcpy(btab + inoffsets[i], &tmp_d, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"h",strlen(members[i])) == 0 ||
                strncmp(members[i],"mfp",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf;
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"vx",strlen(members[i])) == 0 ||
                strncmp(members[i],"vy",strlen(members[i])) == 0 ||
                strncmp(members[i],"vz",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf/timecf;
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"mass",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * masscf;
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"u",strlen(members[i])) == 0 ) {

                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf*lengthcf/(timecf*timecf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"udot",strlen(members[i])) == 0 ) {

                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf*lengthcf/(timecf*timecf*timecf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"rho",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * masscf/(lengthcf*lengthcf*lengthcf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"drho_dt",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * masscf/(lengthcf*lengthcf*lengthcf*timecf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"ax",strlen(members[i])) == 0 ||
                strncmp(members[i],"ay",strlen(members[i])) == 0 ||
                strncmp(members[i],"az",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf/(timecf*timecf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);
                
            } else if( strncmp(members[i],"lax",strlen(members[i])) == 0 ||
                strncmp(members[i],"lay",strlen(members[i])) == 0 ||
                strncmp(members[i],"laz",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf/(timecf*timecf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"phi",strlen(members[i])) == 0 ) {
                
                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf*lengthcf/(timecf*timecf);
                /* I think that's correct. double-check this at some point ~CIE*/
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);

            } else if( strncmp(members[i],"pr",strlen(members[i])) == 0 ) {

                tmp = *(float *)(btab + inoffsets[i]);
                tmp = tmp * lengthcf/(masscf*timecf*timecf);
                memcpy(btab + inoffsets[i], &tmp, SDFtype_sizes[ types[i] ]);
                
            } else {
                continue;
            }
        }

        /*dump the outbtab data into the file now-CIE*/
        fwrite(btab, stride, 1, fp);

    }
    fprintf(frp,"wrote %d abundances to file\n", counter);

    /*and we're done! clean up now -CIE*/
    free(members);
    free(addrs);
    free(strides);
    free(nobjs);
    free(starts);
    free(types);
    free(inoffsets);

    free(btab);
}


