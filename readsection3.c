/*

   PURPOSE:
	to read in large SDF files (10 M particles or more, even though it should
	work for any size ) and shrink them by
	just writing every INCR line but all columns to file.
	This routine reads in the whole file, line-by-line (or rather, INCR lines
	at a time since I couldn't get the single line-by-line to work yet), and
	writes the first of each set of lines (i.e. every INCR line) to file.

   COMPILE:
	with Makefile2. un-comment appropriate lines.
	This routine needs some libraries from the tree code and SDF routines,
	so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
	once for serial use (i.e. without PAROS flag).

   SYNTAX:
	readsection3 <large-file.sdf> <outfile.sdf>

   METHOD:
	@ open <large-file.sdf>
	@ read in INCR lines of all columns with SDFseekreadvecs.
	@ write the first line of the read in data block to <outfile.sdf>

   NOTE:
	assumes that npart from the header contains the actual number of particles
	also, in the past the free'ing of allocated arrays has caused segfaults

*/

#include <stdio.h>
#include <stdlib.h>
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

typedef struct {
    float pos[3];
    float h;
    float mass;
    float rho;
} SPHbody;

int INCR = 25;

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
    writescalars(sdfp, fp);/*writes the header for the scalars (non-structs)*/
    writestructs(sdfp, fp);

    /*fclose(fp);*/
    SDFclose(sdfp);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp)
{
    char input;

    if (argc != 3) {
	fprintf(stderr, "Usage: %s SDFfile outfile \n", argv[0]);
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
    int flag = 0, nflag = 0;
    char **vecs;
    SDF_type_t type;
    datum_t datum;

    nvecs = SDFnvecs(sdfp);/*figure out number of lines in the header, basically-CE*/
    vecs = SDFvecnames(sdfp);/*get the names of the variables/parameters in the header-CE*/

/*go through all of them individually-CE*/
    for (i = 0; i < nvecs; ++i) {
        /*figure out if scalar(=1) or array(!=1)?-CE*/
/*
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs != 1) continue;
*/

        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
        if(flag) continue;

        nflag = 0;
        if (strncmp(vecs[i], "npart", strlen(vecs[i])) == 0) nflag = 1;

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

        /* adjust new particle number */
        if(nflag) datum.i = (int)datum.i/INCR;

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

/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp)
{
    int i, j, nvecs, nmembers;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0, outstride = 0;
    void *btab;
    void **addrs;
    int *inoffsets, *starts, *lines, *strides;
    int flag=0;
    int nlines;
    /*make INCR and nlines user input */

    SDFgetint(sdfp, "npart", &nlines);

/* this does not attempt to load the whole file into memory, does it? -CE */
    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid. But I know that the structure (so far) always has "x" as the first
       member, so I can start counting from there -CE */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        /*get columns of interest*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag=1;
	if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

/*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    lines = (int *)malloc(nmembers * sizeof(int));
    strides = (int *)malloc(nmembers * sizeof(int));

/*one by one, go through the fields in the column, i.e. members of the struct?-CE*/
    flag = 0;
    for (i = 0, stride = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
        if(flag) {
            members[nmembers] = vecs[i];
            lines[nmembers] = 1;//INCR;
            starts[nmembers] = 0;
            types[nmembers] = SDFtype(members[nmembers], sdfp);
            inoffsets[nmembers] = stride;/* offsets (from beginning of 'line'?) of each column
                                            of data (struct member) */
            stride += SDFtype_sizes[types[nmembers]];
            printf("member = %s offset: %d \n",
            members[nmembers],inoffsets[nmembers]);
            ++nmembers;
        }
    }

    /*btab = (void *)malloc(nlines * outstride);*/
    printf("malloc'ing btab: %u\n", INCR*stride);
    btab = (void *)malloc(INCR * stride);

    /*calculate the byte offset in memory to the address of the next member-CE*/
	for (i=0; i<nmembers; i++)
	    addrs[i] = (char *)btab + inoffsets[i];

    printf("reading in %d lines ...\n",nlines);

    /*print the struct declaration part from the header-CE*/
    fprintf(fp, "struct {\n");
    for (i = 0; i < nmembers; ++i) {
        strides[i] = stride;

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
	        fprintf(stderr, "%s: type not supported for member %d\n", members[i], i);
	        exit(-1);
        }
	}
    /*print the struct declaration part from the header-CE*/
    fprintf(fp, "}[%d];\n", nlines/INCR);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    /*try reading in one line at a time */
    for(j = 0; j < nlines; j = j+INCR) {

        for( i = 0; i < nmembers; i++)
            starts[i] = j;

        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, lines, addrs,
		     strides);

        /*dump the outbtab data into the file now-CE*/
        fwrite(btab, stride, 1, fp);

    }


/*and we're done! clean up now -CE: if it ever works*/
    free(members);
    free(btab);
    free(addrs);
    free(types);
    free(inoffsets);
    free(starts);
    free(lines);

}
