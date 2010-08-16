/*

   PURPOSE:
        to merge two SDF files into one file.

   COMPILE:
	with Makefile2. un-comment appropriate lines.
	This routine needs some libraries from the tree code and SDF routines,
	so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
	once for serial use (i.e. without PAROS flag).

   RUN:
	mergeSDFs2 <file-1.sdf> <file-2.sdf> <outfile.sdf>

	You are then prompted for a radius at which the two files should be
        merged, or enter -99.0 to just paste the second file onto the first.

   METHOD:
	The code first checks which of the two headers is shorter (has fewer
	members) and writes the shorter one to <outfile.sdf> (I'm not sure
	why I did it that way, but I assume I had a good reason. Maybe it's
	because it is easier to skip extra data than to read data that is
	not there).
	The code reads in <file-1.sdf> line-by-line, and if a radius of -99.0
	was entered writes the whole file to <oufile.sdf>, else it writes
	all particles with radius less than entered radius to <outfile.sdf>.
	Then <file-2.sdf> is read in ,line-by-line, and added to <outfile.sdf>.
	The code assumes that all particles of <file-2.sdf> have a radius
	greater than the merge radius, so all particles are added.
	The particle ID's of <file-2.sdf> are updated by adding the largest
	ID value from <file-1.sdf>.
	The npart values in the header of <outfile.sdf> are also updated
	to reflect the new total particle number.

   NOTE: the original files are not modified.

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

int GetNumOfDigits(int number);

static void initargs(int argc, char *argv[], SDF **sdfp1, SDF **sdfp2, FILE **fp);
static void writeinit(FILE *fp);
/*static void writescalars(SDF *sdfp1, FILE *fp, fpos_t *pos_npart, int npart);*/
static void writestructs(SDF *sdfp1, SDF *sdfp2, FILE *fp);

int main(int argc, char *argv[])
{
    SDF *sdfp1 = NULL;
    SDF *sdfp2 = NULL;
    FILE *fp = NULL;

    initargs(argc, argv, &sdfp1, &sdfp2, &fp);

    writeinit(fp);
    writestructs(sdfp1, sdfp2, fp);

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
}

static void writescalars(SDF *sdfp, FILE *fp, fpos_t *pos_npart, int npart)
{
    int i, nvecs;
    int flag = 0;
    char **vecs;
    SDF_type_t type;
    datum_t datum;
    fpos_t pos;

    /*figure out number of lines in the header, basically-CE*/
    nvecs = SDFnvecs(sdfp);

    /*get the names of the variables/parameters in the header-CE*/
    vecs = SDFvecnames(sdfp);

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
                datum.i = npart;
            }
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
    *pos_npart = pos;
}

/*this writes the actual data*/
static void writestructs(SDF *sdfp1, SDF *sdfp2, FILE *fp)
{
    int i, j, nvecs1, nvecs2, nmembers, nmembers1, nmembers2, maxmbrs;
    char **vecs1, **vecs2, **members;
    SDF_type_t *types;
    size_t outstride = 0;
    void *btab;
    void **addrs;
    int *inoffsets, *strides, *starts, *lines;
    int incr=1, stride = 0, flag=0;
    int npart1, npart2, newnpart, countnpart;
    int digits, newdigits,ident, idindex, identmax = 0;
    fpos_t pos1_npart, pos_npart;
    double x, y, z;
    float radius, R0;
    int abarr1[2][22], abarr2[2][22];
    int abinfo1, abinfo2, first=1, abundflag=0,k;

    /* count number of data columns, i.e. struct members with more than one element */
    /* this method assumes that "x" is the first data column */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid.
     */

    nvecs1 = SDFnvecs(sdfp1);
    vecs1 = SDFvecnames(sdfp1);

    flag = 0;
    for (i = 0, nmembers1 = 0; i < nvecs1; ++i) {
        if (strncmp(vecs1[i], "x", strlen(vecs1[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
	if (flag) ++nmembers1;
    }

    nvecs2 = SDFnvecs(sdfp2);
    vecs2 = SDFvecnames(sdfp2);

    flag = 0;
    for (i = 0, nmembers2 = 0; i < nvecs2; ++i) {
        if (strncmp(vecs2[i], "x", strlen(vecs2[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
	if (flag) ++nmembers2;
    }

    /* prompt user for merger radius */
    printf("Enter merge radius or -99. : ");
    scanf("%f", &R0);
    printf(" %f\n", R0);

    SDFgetint(sdfp1, "npart", &npart1);
    SDFgetint(sdfp2, "npart", &npart2);
    printf("%d and %d particles\n", npart1, npart2);

    if( nmembers1 <= nmembers2) {
       maxmbrs = nmembers1;
    } else { maxmbrs = nmembers2;}

    /*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(maxmbrs * sizeof(char *));
    addrs = (void **)malloc(maxmbrs * sizeof(void *));
    types = (SDF_type_t *)malloc(maxmbrs * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(maxmbrs * sizeof(int));
    strides = (int *)malloc(maxmbrs * sizeof(int));
    starts = (int *)malloc(maxmbrs * sizeof(int));
    lines = (int *)malloc(maxmbrs * sizeof(int));

    printf("done malloc'ing\n");

    /* find the file with the fewer members and write its header to output */
    /* done by populating the output arrays from the appropriate file */

    flag = 0;
    if ( nmembers1 < nmembers2) {
    /*one by one, go through the fields in the column, i.e. members of the struct?-CE*/
        for (i = 0, stride = 0, nmembers = 0; i < nvecs1; ++i) {

            if (strncmp(vecs1[i], "x", strlen(vecs1[i])) == 0) flag=1;

            if(flag) {
	        members[nmembers] = vecs1[i];
	        types[nmembers] = SDFtype(members[nmembers], sdfp1);
	        inoffsets[nmembers] = stride;
	        stride += SDFtype_sizes[types[nmembers]];
                lines[nmembers] = incr;
                nmembers++;
            }

            /* get index that holds 'ident' for later updating */
            if (strncmp(vecs1[i], "ident", strlen(vecs1[i])) == 0) idindex=nmembers-1;
            if(strncmp(vecs1[i],"p1",strlen(vecs1[i])) == 0) abinfo1=nmembers-1;

        }
        for (i = 0; i < nmembers; i++)
            strides[i] = stride;

    } else {
    /*one by one, go through the fields in the column, i.e. members of the struct?-CE*/
        for (i = 0, stride = 0, nmembers = 0; i < nvecs2; ++i) {

            if (strncmp(vecs2[i], "x", strlen(vecs2[i])) == 0) flag=1;

            if(flag) {
	        members[nmembers] = vecs2[i];
	        types[nmembers] = SDFtype(members[nmembers], sdfp2);
	        inoffsets[nmembers] = stride;
	        stride += SDFtype_sizes[types[nmembers]];
                lines[nmembers] = incr;
                nmembers++;
	    }

            /* get index that holds 'ident' for later updating */
            if (strncmp(vecs1[i], "ident", strlen(vecs1[i])) == 0) idindex=nmembers-1;
            if(strncmp(vecs2[i],"p1",strlen(vecs2[i])) == 0) abinfo2=nmembers-1;

        }
        for (i = 0; i < nmembers; i++)
            strides[i] = stride;

    }

    /* unnecesary, just use 'stride' ? CE */
    outstride = 0;
    for(i=0; i< maxmbrs; i++) outstride += SDFtype_sizes[ types[i] ];

    printf("outstride = %d\n", outstride);

    /* holds the read in data */
    btab = (void *)malloc( stride*incr );


    printf("printing header\n");

    /* writes the new header from the appropriate file. Assumes that both files
     * are merged in their entirety. Returns the location in memory where the
     * total particle number is stored so it can be updated later */

    if (nmembers1 < nmembers2) {
        writescalars(sdfp2, fp, &pos_npart,npart1+npart2);
    } else {
        writescalars(sdfp1, fp, &pos_npart,npart1+npart2);
    }

    /*print the struct declaration part from the header-CE*/
    fprintf(fp, "struct {\n");
    for (i = 0; i < maxmbrs; ++i) {
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
    /* again, assume that all of both files is merged. Save this memory location also */
    fgetpos(fp, &pos1_npart);
    fprintf(fp, "}[%d];\n", npart1+npart2);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    /*calculate the byte offset in memory to the address of the next member-CE*/
    addrs[0] = (char *)btab;
    for (i=1; i< maxmbrs; i++)
         addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];

    /*loop over each file consecutively in chunks of data, write that chunk to file*/
    printf("getting file 1 .... ");
    for( i=0, countnpart = 0; i < npart1; i++) {
        /* need to increment starts-array */
        for( j = 0; j < maxmbrs; j++)
            starts[j] = i;

        /* read data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp1, maxmbrs, members, starts, lines, addrs, strides);

        ident = *((int *)(btab + inoffsets[idindex]));

        /* calculate radius for merging */
        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        radius = sqrt(x*x + y*y + z*z);

        if(first) {
            for(k=0;k<22;k++){
                abarr1[0][k] = *((int *)(btab + inoffsets[abinfo1+k]));/*Z*/
                abarr1[1][k] = *((int *)(btab + inoffsets[abinfo1+22+k]));/*N*/
            }
        }
        first=0;

	/*dump the btab data into the file now-CE*/
        if( (radius < R0 ) ) {
	    fwrite(btab, outstride, 1, fp);
            /* count number of particles actually writted */
            countnpart++;
            if(ident > identmax)
               identmax = ident;
        }
        else if( R0 < 0.0 ) {
            fwrite(btab, outstride, 1, fp);
            /* count number of particles actually writted */
            countnpart++;
            if(ident > identmax)
               identmax = ident;
        }

    }
    newnpart = countnpart;

    first=1;

    printf("got %d lines\n",countnpart);
    printf("getting file 2 .... ");
    for( i=0, countnpart = 0; i < npart2; i++) {
        /* need to increment starts-array */
        for( j = 0; j < maxmbrs; j++)
            starts[j] = i;

        /* read data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp2, maxmbrs, members, starts, lines, addrs, strides);

        /* update the particle id, so it does not start at 1 again */
        ident = *((int *)(btab + inoffsets[idindex]));
        ident += identmax + 1;
        memcpy( btab + inoffsets[idindex], &ident, sizeof(ident) );

        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        radius = sqrt(x*x + y*y + z*z);

        /* check if the abundance information is the same in both files, and spit
         * out a warning if it is not */
        /* Note: doesn't seem to be working yet */
        if(first) {
            for(k=0;k<22;k++){
                abarr2[0][k] = *((int *)(btab + inoffsets[abinfo2+k]));/*Z*/
                abarr2[1][k] = *((int *)(btab + inoffsets[abinfo2+22+k]));/*N*/
                if((abarr1[0][k] != abarr2[0][k]) || (abarr1[1][k] != abarr2[1][k]))
                   abundflag=1;
            }
        }
        if(first && (abundflag == 1)) {
           printf("warning: abundance info in files might not be the same!\n");
           printf("%10s %10s\n%5s%5s %5s%5s\n", "file1", "file2","p","n","p","n");
           for(k=0;k<22;k++) printf("%5d%5d %5d%5d\n",
               abarr1[0][k],abarr1[1][k],abarr2[0][k],abarr2[1][k]);
        }
        first=0;
	/*dump the btab data into the file now-CE*/
        if( radius > R0 ) { /* radius is always greater than -99., no extra case needed*/
	    fwrite(btab, outstride, 1, fp);
            /* count number of particles actually writted */
            countnpart++;
        }

    }
    newnpart += countnpart;
    printf("got %d lines\n",countnpart);

    newdigits = GetNumOfDigits(newnpart);
    digits = GetNumOfDigits(npart1+npart2);

    /* update npart to the new value */
    fsetpos(fp, &pos1_npart);
    fprintf(fp, "}[%d];", newnpart);

    /*add blanks in case newnpart has fewer digits than npart1+npart2*/
    while(digits > newdigits) {
        fprintf(fp, " ");
        digits--;
    }

    digits = GetNumOfDigits(npart1+npart2);

    fsetpos(fp, &pos_npart);
    fprintf(fp, "int %s = %d;", "npart", newnpart);
    while(digits > newdigits){
        fprintf(fp, " ");
        digits--;
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


int GetNumOfDigits(int number) {
    int digits = 1;

    while(number/10 >= 1) {
        number = number/10;
        digits++;
    }
    return digits;
}
