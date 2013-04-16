/*
   PURPOSE:
    to add a perturbation to some particles to an sdf file
        This version should also work for large files, since it reads in the
        sdf file line-by-line.

   COMPILE:
    with Makefile.
    make ARCH=<arch> CC=<c-compilter> PROGS=addperturb -f Makefile
    This routine needs some libraries from the tree code and SDF routines,
    so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
    once for serial use (i.e. without PAROS flag).

   RUN:
    addvel <in-file.sdf> <out-file.sdf>

   METHOD:
    @ read in header of <in-file.sdf>
    @ prompt user for the shock radius (or below where-ever the velocity
      asymmetry is to be added
    @ read in data line-by-line
    @ calculate radius for each particle
    @ determine if particle is in/inside of the shock (by simple comparison
      to the input radius)
    @ if yes, modify velocity according to prescription
    @ write data to <out-file.sdf>

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

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp);
static void writeinit(FILE *fp);
static void writescalars(SDF *sdfp, FILE *fp);
static void writestructs(SDF *sdfp, FILE *fp);

int perturb_which = 0;


double get_ye(float abundarr[], int nparr[], int nnarr[], int Niso);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;

    initargs(argc, argv, &sdfp, &fp);

    writeinit(fp);
    writescalars(sdfp, fp);
    writestructs(sdfp, fp);

    fclose(fp);
    SDFclose(sdfp);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp)
{
    char input;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s SDFfile outfile perturb-switch\n", argv[0]);
        exit(1);
    }

    *sdfp = SDFopen(NULL, argv[1]);
    if (*sdfp == NULL) {
        fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
        exit(2);
    }

    if (access(argv[2], F_OK) == 0) {
        fprintf(stderr, "%s: %s exists; overwrite (y/n)? ", argv[0],
            argv[2]);
        input = getc(stdin);
        if ((input != 'y') && (input != 'Y'))
            exit(3);
    }
    *fp = fopen(argv[2], "w");
    if (*fp == NULL) {
        fprintf(stderr, "%s: %s\n", argv[2], strerror(errno));
        exit(errno);
    }
   
    perturb_which = atoi(argv[3]);
    if(perturb_which == 0) printf("adding random perturbation to density (via h)\n");
    else if(perturb_which == 1) printf("adding random perturbation to velocity\n");
    else if(perturb_which == 2) printf("adding random perturbation to density (via h) and velocity\n");
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
    int ixvel, ih;
    int *strides, *nobjs, *starts, *inoffsets;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0;
    void *btab;
    void **addrs;
    double alpha, beta, vfactor;
    float set_radius;
    double x, y, z, radius;
    float vx, vy, vz, vx2, vy2, vz2, h, h3, perturb;

    frp = fopen("log.out", "w");
    if(!frp) printf("error opening log file!\n");

    printf("set_radius: enter 0 if set by vr:");
    scanf("%f",&set_radius);
    printf("%f\n",set_radius);

    vfactor = 1.00;

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
            if (strncmp(vecs[i], "vx", strlen(vecs[i])) == 0) ixvel = nmembers;
            if (strncmp(vecs[i], "h", strlen(vecs[i])) == 0) ih = nmembers;
            stride += SDFtype_sizes[types[nmembers]];

            ++nmembers;
        }
    }

    fprintf(frp,"vx at %d\n",ixvel);
    fprintf(frp,"h at %d\n",ih);

    btab = (void *)malloc(stride);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
    for (i = 0; i < nmembers; ++i) {
        addrs[i] = (char *)btab + inoffsets[i];
        strides[i] = stride;
    }


/* up to here just the SDF file is read in, so everything stays the same -CIE */

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
        //printf("%8d ",j);
        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
                     strides);

        /* calculate radius. this should work -CIE */
        x = *(double *)(btab + inoffsets[0]);
        y = *(double *)(btab + inoffsets[1]);
        z = *(double *)(btab + inoffsets[2]);
        radius = sqrt( x*x + y*y + z*z );

        /* get velocities and density */
        vx = *(float *)(btab + inoffsets[ixvel]);
        vy = *(float *)(btab + inoffsets[ixvel+1]);
        vz = *(float *)(btab + inoffsets[ixvel+2]);
        h = *(float *)(btab + inoffsets[ih]);

        /* calculate the velocity asymmetry */
        /* only for expanding velocities */
        if( ((vx*x+vy*y+vz*z)/radius > 0.0 && set_radius == 0.)
            || (radius < set_radius)) {

            if(perturb_which==0 || perturb_which == 2){ /* perturbing velocity */
                perturb = (double)rand()/(double)RAND_MAX-0.5;
                vx2 = vx*(1.+perturb/5.);
                vy2 = vy*(1.+perturb/5.);
                vz2 = vz*(1.+perturb/5.);
                memcpy(btab + inoffsets[ixvel], &vx2, SDFtype_sizes[ types[ixvel] ]);
                memcpy(btab + inoffsets[ixvel+1], &vy2, SDFtype_sizes[ types[ixvel+1] ]);
                memcpy(btab + inoffsets[ixvel+2], &vz2, SDFtype_sizes[ types[ixvel+2] ]);

            } 

            if(perturb_which == 1 || perturb_which == 2) { /* perturbing density */

                perturb = (double)rand()/(double)RAND_MAX-0.5;
                h3 = h*h*h;
                h3 = h3* (1.+perturb/5.);
                h = pow(h, 0.33333333333333333333);
                memcpy(btab + inoffsets[ih ], &h, SDFtype_sizes[ types[ ih ] ]);

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


