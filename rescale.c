/*
   PURPOSE:
	to add a velocity to an sdf file
        This version should also work for large files, since it reads in the
        sdf file line-by-line.

   COMPILE:
    with Makefile2. un-comment appropriate lines.
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

    if (argc != 3) {
	fprintf(stderr, "Usage: %s SDFfile outfile\n", argv[0]);
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
/* 	    exit(-1); */
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
/* 	    exit(-1); */
	}
    }
}

static void writestructs(SDF *sdfp, FILE *fp)
{
    int i, j, nvecs, nmembers;
    int len, counter, flag = 0;
    int nrecs;
    int ixvel;
    int *strides, *nobjs, *starts, *inoffsets;
    int mflag=0, hindex, rhoindex, vflag=0, aflag=0;
    char **vecs, **members;
    char member[10];
    SDF_type_t *types;
    size_t stride = 0;
    void *btab;
    void **addrs;
    float set_radius, factor=4.5;
    double x, y, z, radius;
    float var, var_new, var_ave, rho, h;


    printf("quantity to rescale: \n");
    scanf("%s",member);
    if( strncmp(member, "mass", 4) == 0) mflag = 1;
    if( strncmp(member, "vel", 3) == 0) vflag = 1;
    if( strncmp(member, "acc", 3) == 0) aflag = 1;

    printf("resetting %s above r=\n", member);
    scanf("%f",&set_radius);

    printf("by factor:\n");
    scanf("%f",&factor);

    printf("around average value:\n");
    scanf("%f", &var_ave);

    SDFgetint(sdfp, "npart", &nrecs);
    printf("%d particles\n", nrecs);

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) ++nmembers;
    }


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
            if (strncmp(vecs[i], member, strlen(vecs[i])) == 0) ixvel = nmembers;
	    stride += SDFtype_sizes[types[nmembers]];

        if(mflag) {
            if (strncmp(vecs[i], "h", strlen(vecs[i])) == 0) hindex = nmembers;
            if (strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) rhoindex = nmembers;
        }
        if( vflag ) {
            if( strncmp(vecs[i], "vx", strlen(vecs[i])) == 0) ixvel = nmembers;
        }
        if( aflag ) {
            if( strncmp(vecs[i], "ax", strlen(vecs[i])) == 0) ixvel = nmembers;
        }
	    ++nmembers;
	}
    }

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
        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
		     strides);

        /* calculate radius. this should work -CIE */
        x = *(double *)(btab + inoffsets[0]);
        y = *(double *)(btab + inoffsets[1]);
        z = *(double *)(btab + inoffsets[2]);
        radius =  sqrt(x*x + y*y + z*z);
        if(mflag)
        rho = *(float *)(btab + inoffsets[rhoindex]);

        if( radius >= set_radius) {

            /* get quantity */
            var = *(float *)(btab + inoffsets[ixvel]);
            var_new = (var-var_ave)*factor+var_ave;
    
            if(mflag) {
                h = *(float *)(btab + inoffsets[hindex]);
                //h = pow( (var_new/(0.6214123964*rho)), 0.333333333333333333);
                h = h*pow( (var_new/var), 0.33333333333333 );
            }

            memcpy(btab + inoffsets[ixvel], &var_new, SDFtype_sizes[ types[ixvel] ]);

            if(mflag)
                memcpy(btab + inoffsets[hindex], &h, SDFtype_sizes[ types[hindex] ]);

            if( aflag ) {
                /* get quantity - ay*/
                var = *(float *)(btab + inoffsets[ixvel+1]);
                var_new = (var-var_ave)*factor+var_ave;
                memcpy(btab + inoffsets[ixvel+1], &var_new, SDFtype_sizes[ types[ixvel] ]);
    
                /* get quantity - az */
                var = *(float *)(btab + inoffsets[ixvel+2]);
                var_new = (var-var_ave)*factor+var_ave;
                memcpy(btab + inoffsets[ixvel+2], &var_new, SDFtype_sizes[ types[ixvel] ]);
            }

            if( vflag ) {
                /* get quantity - vy*/
                var = *(float *)(btab + inoffsets[ixvel+1]);
                var_new = (var-var_ave)*factor+var_ave;
                memcpy(btab + inoffsets[ixvel+1], &var_new, SDFtype_sizes[ types[ixvel] ]);
    
                /* get quantity - vz */
                var = *(float *)(btab + inoffsets[ixvel+2]);
                var_new = (var-var_ave)*factor+var_ave;
                memcpy(btab + inoffsets[ixvel+2], &var_new, SDFtype_sizes[ types[ixvel] ]);
            }
        }

        /*dump the outbtab data into the file now-CIE*/
        fwrite(btab, stride, 1, fp);

    }
    printf("wrote %d abundances to file\n", counter);

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
