/*
   PURPOSE:
        to merge two SDF files into one file.
        You can now enter the radius at which the two files should be
        merged, or -99. to just paste the second file onto the first.

   NOTE: the original files are not modified.
         - assumes currently that the body is the same, and the header from the 
           first file is written

   DONE: 1) get list of columns in both SDF files
   DONE: 2) compare and exit if the headers are not the same
   DONE: 2) read all scalars
   DONE: 3) read all structure members
   DONE: 5) loop over 3-4 until whole file is read (or seg fault is reached ;-P)
   DONE: 6) write scalars and buffer to file
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <SDF.h>

typedef enum SDF_type_enum SDF_type_t;

typedef union {
    int i;
    float f;
    double d;
} datum_t;

static void initargs(int argc, char *argv[], SDF **sdfp1, SDF **sdfp2, FILE **fp);
static void writeinit(FILE *fp);
static void writescalars(SDF *sdfp1, FILE *fp, fpos_t *pos_npart, int max_npart);
static void writestructs(SDF *sdfp1, SDF *sdfp2, FILE *fp, fpos_t pos_npart);

int main(int argc, char *argv[])
{
    SDF *sdfp1 = NULL;
    SDF *sdfp2 = NULL;
    FILE *fp = NULL;
    fpos_t pos_npart;
    int max_npart, npart1, npart2;

    initargs(argc, argv, &sdfp1, &sdfp2, &fp);

    writeinit(fp);

    /* determine the max number of particles */
    SDFgetint(sdfp1, "npart", &npart1);
    SDFgetint(sdfp2, "npart", &npart2);
    max_npart = npart1+npart2;

    writescalars(sdfp1, fp, &pos_npart, max_npart);/*writes the header for the scalars (non-structs)*/
    writestructs(sdfp1, sdfp2, fp, pos_npart);

    fclose(fp);
    SDFclose(sdfp1);
    SDFclose(sdfp2);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp1, SDF **sdfp2, FILE **fp)
{
    char input;

    if (argc != 4) {
	fprintf(stderr, "Usage: %s SDFfile SDFfile outfile \n", argv[0]);
	exit(1);
    }

    *sdfp1 = SDFopen(NULL, argv[1]);
    if (*sdfp1 == NULL) {
	fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
	exit(2);
    }

    *sdfp2 = SDFopen(NULL, argv[2]);
    if (*sdfp2 == NULL) {
	fprintf(stderr, "%s: %s: %s\n", argv[0], argv[2], SDFerrstring);
	exit(2);
    }

    if (access(argv[3], F_OK) == 0) {
        fprintf(stderr, "%s: %s exists; overwrite (y/n)? ", argv[0],
		argv[3]);
        input = getc(stdin);
        if ((input != 'y') && (input != 'Y'))
            exit(3);
    }

    *fp = fopen(argv[3], "w");
    if (*fp == NULL) {
	fprintf(stderr, "%s: %s\n", argv[3], strerror(errno));
	exit(errno);
    }

}

static void writeinit(FILE *fp)
{
    fpos_t pos1;
    fprintf(fp, "# SDF\n");
    fprintf(fp, "parameter byteorder = %#x;\n", SDFcpubyteorder());
    fgetpos(fp, &pos1);
    printf("hello\n");
}

static void writescalars(SDF *sdfp, FILE *fp, fpos_t *pos_npart, int max_npart)
{
    int i, nvecs;
    int flag = 0;
    char **vecs;
    SDF_type_t type;
    datum_t datum;
    fpos_t pos;

    nvecs = SDFnvecs(sdfp);/*figure out number of lines in the header, basically-CE*/
    vecs = SDFvecnames(sdfp);/*get the names of the variables/parameters in the header-CE*/

/*go through all of them individually-CE*/
    for (i = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) continue;

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
            if( !strncmp(vecs[i], "npart", strlen(vecs[i])) ) {
                fgetpos(fp, &pos);
	        fprintf(fp, "int %s = %d;\n", vecs[i], max_npart);
            } else {
	        fprintf(fp, "int %s = %d;\n", vecs[i], datum.i);
            }
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
    *pos_npart = pos;
}

/*this writes the actual data*/
static void writestructs(SDF *sdfp1, SDF *sdfp2, FILE *fp, fpos_t pos_npart)
{
    int i, j, nvecs1, nvecs2, nmembers1, nmembers2, nmembers;
    char **vecs1, **vecs2, **members1, **members2;
    SDF_type_t *types;
    SDF_type_t type2;
    size_t stride = 0, outstride = 0;
    void *btab;
    void **addrs;
    int *inoffsets, *strides, *starts, *lines;
    int flag=0, incr=1;
    int npart1, npart2, newnpart, countnpart, max;
    int ident, idindex, ident_last, ident_max, identnew;
    fpos_t pos1_npart;
    double x, y, z;
    float radius, R0, maxR0;
    /*make INCR and nlines user input */

    nvecs1 = SDFnvecs(sdfp1);
    vecs1 = SDFvecnames(sdfp1);
    max = nvecs1;

    nvecs2 = SDFnvecs(sdfp2);
    vecs2 = SDFvecnames(sdfp2);

    /* Count structure members */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid. But I know that the structure (so far) always has "x" as the first
       member, so I can start counting from there -CE */
    for (i = 0, nmembers1 = 0; i < nvecs1; ++i) {
        if (strncmp(vecs1[i], "x", strlen(vecs1[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
	if (flag) ++nmembers1;
    }

    flag=0;
    for (i = 0, nmembers2 = 0; i < nvecs2; ++i) {
        if (strncmp(vecs2[i], "x", strlen(vecs2[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
	if (flag) ++nmembers2;
    }

    if( nmembers1 != nmembers2) {
        fprintf(stderr, "non-matching headers: %d and %d\n",nvecs1,nvecs2);
//        exit(1);
    }

    printf("Enter merge radius or -99. : ");
    scanf("%f", &R0);
    printf(" %f\n", R0);

    printf("Enter max r: ");
    scanf("%f",&maxR0);
    printf("%f\n",maxR0); 

    SDFgetint(sdfp1, "npart", &npart1);
    SDFgetint(sdfp2, "npart", &npart2);
    printf("%d and %d particles\n", npart1, npart2);

/*malloc memory space for the respective features of the struct-CE*/
    members1 = (char **)malloc(nmembers1 * sizeof(char *));
    members2 = (char **)malloc(nmembers2 * sizeof(char *));
    addrs = (void **)malloc(nmembers1 * sizeof(void *));
    types = (SDF_type_t *)malloc(nmembers1 * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers1 * sizeof(int));
    strides = (int *)malloc(nmembers1 * sizeof(int));
    starts = (int *)malloc(nmembers1 * sizeof(int));
    lines = (int *)malloc(nmembers1 * sizeof(int));

    printf("done malloc'ing\n");
    flag=0;
/*one by one, go through the fields in the column, i.e. members of the struct?-CE*/
    for (i = 0, stride = 0, nmembers = 0; i < nvecs1; ++i) {
        if (strncmp(vecs1[i], "x", strlen(vecs1[i])) == 0) flag=1;
        if (strncmp(vecs1[i], "ident", strlen(vecs1[i])) == 0) idindex=nmembers;
        if(flag) {
	    members1[nmembers] = vecs1[i];
//	    members2[nmembers] = vecs2[i];
	    types[nmembers] = SDFtype(members1[nmembers], sdfp1);
	    inoffsets[nmembers] = stride;/* offsets (from beginning of 'line'?) of each column
                                            of data (struct member) */
	    stride += SDFtype_sizes[types[nmembers]];
            lines[nmembers] = incr;
            nmembers++;
	}
    }

    /* unnecesary, just use 'stride' ? CE */
    outstride = 0;
    for(i=0; i< nmembers1; i++) outstride += SDFtype_sizes[ types[i] ];

    printf("outstride = %d\n", outstride);

    btab = (void *)malloc( outstride*incr );

    printf("printing header\n");
    /*print the struct declaration part from the header-CE*/
    fprintf(fp, "struct {\n");
    for (i = 0; i < nmembers1; ++i) {
        strides[i] = stride;
	switch (types[i]) {
	case SDF_INT:
	    fprintf(fp, "\tint %s;\n", members1[i]);
	    break;
	case SDF_FLOAT:
            fprintf(fp, "\tfloat %s;\n", members1[i]);
	    break;
	case SDF_DOUBLE:
	    fprintf(fp, "\tdouble %s;\n", members1[i]);
	    break;
	default:
	    fprintf(stderr, "%s: type not supported\n", members1[i]);
	    exit(-1);
	}
    }
    fgetpos(fp, &pos1_npart);
    fprintf(fp, "}[%d];\n", npart1+npart2); /*figure out how to save the location of the file position
                                       indicator, so we can come back and update npart -CE */
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    /*calculate the byte offset in memory to the address of the next member-CE*/
    addrs[0] = (char *)btab;
    for (i=1; i< nmembers1; i++) addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];

/*loop over each file consecutively in chunks of data, write that chunk to file*/

    printf("getting file 1 .... ");
    for( i=0, countnpart = 0; i < npart1-1; i++) {
        /* need to increment starts-array */
        for( j = 0; j < nmembers1; j++)
            starts[j] = i;

        /* read data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp1, nmembers1, members1, starts, lines, addrs, strides);

        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        radius = sqrt(x*x + y*y + z*z);

	/*dump the btab data into the file now-CE*/
        if( radius < R0 ) {
            ident_last=ident;
            ident = *((int *)(btab + inoffsets[idindex]));
            if(ident_last==ident)printf("same ident\n");
	    fwrite(btab, outstride, 1, fp);
            countnpart++;
        }
        else if( R0 == -99.0 ) {
            ident_last=ident;
            ident = *((int *)(btab + inoffsets[idindex]));
            if(ident_last==ident)printf("same ident\n");
            fwrite(btab, outstride, 1, fp);
            countnpart++;
        }

    }
    newnpart = countnpart-1;
    ident_max=ident;

    printf("got %d lines, last ident=%d\n",countnpart,ident_max);
    printf("getting file 2 .... ");
    for( i=0, countnpart = 0; i < npart2-1; i++) {
        /* need to increment starts-array */
        for( j = 0; j < nmembers1; j++)
            starts[j] = i;

        /* read data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp2, nmembers1, members1, starts, lines, addrs, strides);

        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        radius = sqrt(x*x + y*y + z*z);

	/*dump the btab data into the file now-CE*/
        if( radius > R0 && radius <= maxR0) { /* radius is always greater than -99., no extra case needed*/
            /* update the particle id, so it does not start at 1 again */
            ident_last=ident;
            ident = *((int *)(btab + inoffsets[idindex]));
            if(ident_last==ident)printf(" same ident ");
            ++countnpart;
            identnew = countnpart + ident_max; /* so there aren't duplicate particle ids */
            memcpy( btab + inoffsets[idindex], &identnew, sizeof(ident) );

	    fwrite(btab, outstride, 1, fp);
        } 

    }
    newnpart += countnpart-1;
    printf("got %d lines\n",countnpart);

    /* update npart to the new value, adjust for differing number of digits */
    if( (npart1+npart2)/newnpart >= 100 ) {
        fsetpos(fp, &pos1_npart);
        fprintf(fp, "} [%d]; ", newnpart);
    
        fsetpos(fp, &pos_npart);
        fprintf(fp, "int %s = %d;  ", "npart", newnpart);
    } else if( (npart1+npart2)/newnpart >= 10 ) {
        fsetpos(fp, &pos1_npart);
        fprintf(fp, "}[%d]; ", newnpart);
    
        fsetpos(fp, &pos_npart);
        fprintf(fp, "int %s = %d; ", "npart", newnpart);
    } else {
        fsetpos(fp, &pos1_npart);
        fprintf(fp, "}[%d];", newnpart);
    
        fsetpos(fp, &pos_npart);
        fprintf(fp, "int %s = %d;", "npart", newnpart);
    }

/*and we're done! clean up now -CE: if it ever works*/
/*    free(members);
    free(btab);
    free(addrs);
    free(types);
    free(inoffsets);
*/
    /*free(outbtab);*/
}
