/*
   todoublepos: convert positions in SDF files from float to double;
       copy other file contents to output unchanged

   DONE: 1) get list of names in SDF file
   DONE: 2) read all scalars
   DONE: 3) read all structure members, convert float to doubles as appropriate
   DONE: 4) calculate offsets, write adjusted members into buffer
   DONE: 5) write scalars and buffer to file
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
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
    char **vecs;
    SDF_type_t type;
    datum_t datum;

    nvecs = SDFnvecs(sdfp);/*figure out number of lines in the header, basically-CE*/
    vecs = SDFvecnames(sdfp);/*get the names of the variables/parameters in the header-CE*/

/*go through all of them individually-CE*/
    for (i = 0; i < nvecs; ++i) {
	nrecs = SDFnrecs(vecs[i], sdfp);/*figure out if scalar(=1) or array(!=1)?-CE*/
	if (nrecs != 1) continue;

	type = SDFtype(vecs[i], sdfp);
	/* Doesn't handle string */
	/* Also, SDFgetfloat was happy to read a double scalar and
	   convert it for me; that probably works backwards too.  I
	   don't think there's an equivalent for the SDFrdvecs family
	   though, so "read; convert type; write" is the general
	   path. */
/*read in header file, line by line, with the appropriate function-CE*/
	switch (type) {
	case SDF_INT:
/*SDFget*:sdfp=sdf file; vecs[i]= variable name; datum.*=holds value of that variable-CE*/
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

/*write header file, line by line, as the appropriate data type-CE*/
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
    int i, j, nvecs, nmembers, nrecs;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0, outstride = 0;
    void *btab, *outbtab;
    void **addrs;
    int *strides, *nobjs, *starts, *inoffsets, *outoffsets;
    double x, y, z;
    float radius, mass;

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs > 1) ++nmembers;
    }

/*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    strides = (int *)malloc(nmembers * sizeof(int));
    nobjs = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));
    outoffsets = (int *)malloc(nmembers * sizeof(int));

/*one by one, go through all the fields in the column, i.e. members of the struct?-CE*/
    for (i = 0, nmembers = 0, stride = 0; i < nvecs; ++i) {
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs > 1) {
	    members[nmembers] = vecs[i];
	    nobjs[nmembers] = nrecs;/*nobjs[0] is the number of particles. are all elements in
                                      this nobjs array the same then??-CE:yes, but can be
                                      different*/
	    starts[nmembers] = 0;  /* Not correct in parallel; use
				      NobjInitial */
	    types[nmembers] = SDFtype(members[nmembers], sdfp);
	    inoffsets[nmembers] = stride;
	    stride += SDFtype_sizes[types[nmembers]];

	    ++nmembers;
	}
    }

    btab = (void *)malloc(nobjs[0] * stride);

/*calculate the byte offset in memory to the address of the next member-CE*/
    for (i = 0; i < nmembers; ++i) {
	addrs[i] = (char *)btab + inoffsets[i];
	strides[i] = stride;
    }

/*what does this do????-CE*/
    SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
		     strides);

    for (i = 0, outstride = 0; i < nmembers; ++i) {
/*calculate at what byte-intervals the data should be written-CE*/
	outoffsets[i] = outstride;
	outstride += SDFtype_sizes[ types[i] ];
    }

/*malloc enough space in memory for the array that holds the whole output data-CE*/
    outbtab = (void *)malloc(nobjs[0] * outstride);

    for (j = 0; j < nobjs[0]; ++j) {
        x = *((double *)(btab + j*stride + inoffsets[0]));
        y = *((double *)(btab + j*stride + inoffsets[1]));
        z = *((double *)(btab + j*stride + inoffsets[2]));
        radius = sqrt( x*x + y*y + z*z );

	for (i = 0; i < nmembers; ++i) {
            if( !strncmp(members[i], "mass", strlen( members[i] )) ) {
	        mass = *((float *)(btab + j*stride + inoffsets[i]));

                if( (radius < 3.1e-2) && (mass > 5.0) ) {
                    printf("resizing mass %G at r= %G ", mass, radius);
                    mass = mass/2.;
                    printf("to %G\n", mass);
		    memcpy(outbtab + j*outstride + outoffsets[i],
		           &mass, SDFtype_sizes[types[i]]);
                }
                else {
		    memcpy(outbtab + j*outstride + outoffsets[i],
		           btab + j*stride + inoffsets[i],
		           SDFtype_sizes[types[i]]);
                }

            }
            else {
		memcpy(outbtab + j*outstride + outoffsets[i],
		       btab + j*stride + inoffsets[i],
		       SDFtype_sizes[types[i]]);
            }
	}
    }

/*print the struct declaration part from the header-CE*/
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
    fprintf(fp, "}[%d];\n", nobjs[0]);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

/*dump the outbtab data into the file now-CE*/
    fwrite(outbtab, outstride, nobjs[0], fp);

/*and we're done! clean up now -CE*/
    free(members);
    free(addrs);
    free(strides);
    free(nobjs);
    free(starts);
    free(types);
    free(inoffsets);
    free(outoffsets);

    free(btab);
    free(outbtab);
}
