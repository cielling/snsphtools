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

int prep_SDF_args(SDF *sdfp, char **vecs, int nvecs, char ***members, int *nmembers, SDF_type_t **types, int **offsets, int **strides, int **lines, size_t *stride);

int verify_head(char **vecs1, int nvecs1, int nmembers1, char **vecs2, int nvecs2, int nmembers2, int **ordhead);

int verify_body(char **vecs1, int nvecs1, int nmembers1, char **vecs2, int nvecs2, int nmembers2, int **ordbody);

int verify_abun(int abarr1[][22], int abarr2[][22], int nmembers1);




int main(int argc, char *argv[])
{
    SDF *sdfp1 = NULL;
    SDF *sdfp2 = NULL;
    FILE *fp = NULL;

    initargs(argc, argv, &sdfp1, &sdfp2, &fp);

    writeinit(fp);
    writestructs(sdfp1, sdfp2, fp);

    fclose(fp);
/*
    SDFclose(sdfp1);
    SDFclose(sdfp2);
*/

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
    int i, j, k, nvecs1, nvecs2;
	int nmembers, nmembers1, nmembers2, maxmbrs, minmbrs;
    char **vecs1, **vecs2, **members;
    SDF_type_t *types, *outtypes;
	SDF *sdfp = NULL;
    size_t outstride = 0, stride = 0;
    void *btab, *outbtab;
    void **addrs;
    int *inoffsets, *strides, *starts, *lines, *ordhead, *ordbody;
	int *outoffsets;
    int incr=1, flag=0, ind = 0, second = 0;
    int npart, npart1, npart2, newnpart, countnpart;
    int digits, newdigits,ident, idindex, identmax = 0;
	int sort_head=0, sort_body=0;
    fpos_t pos1_npart, pos_npart;
    double x, y, z;
    float radius, R0;
    int abarr1[2][22], abarr2[2][22];
    int abinfo1= -1, abinfo2=-1;

    /* count number of data columns, i.e. struct members with more than one element */
    /* this method assumes that "x" is the first data column */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid.
     */

    nvecs1 = SDFnvecs(sdfp1);
    vecs1 = SDFvecnames(sdfp1);

	ind = 0;

    flag = 0;
    for (i = 0, nmembers1 = 0; i < nvecs1; ++i) {
		//printf("%s\n",vecs1[i]);
        if (strncmp(vecs1[i], "x", strlen(vecs1[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
		if (flag) ++nmembers1;
        if( strncmp(vecs1[i],"p1",strlen(vecs1[i])) == 0) abinfo1=i;
		if( (abinfo1 > 0) && (i-abinfo1 < 44) ) {
			SDFgetint(sdfp1, vecs1[i], &abarr1[ind][(int)((i-abinfo1)/2)]);
			ind = ind > 0 ? 0 : 1;
		}
    }

    nvecs2 = SDFnvecs(sdfp2);
    vecs2 = SDFvecnames(sdfp2);

	ind = 0;

    flag = 0;
    for (i = 0, nmembers2 = 0; i < nvecs2; ++i) {
        if (strncmp(vecs2[i], "x", strlen(vecs2[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
	if (flag) ++nmembers2;
        if( strncmp(vecs2[i],"p1",strlen(vecs2[i])) == 0) abinfo2=i;
		if( (abinfo2 > 0) && (i-abinfo2 < 44) ) {
			SDFgetint(sdfp2, vecs2[i], &abarr2[ind][(int)((i-abinfo2)/2)]);
			ind = (ind != 1);
		}
    }

    /* prompt user for merger radius */
    printf("Enter merge radius or -99. : ");
    scanf("%f", &R0);
    printf(" %f\n", R0);

    SDFgetint(sdfp1, "npart", &npart1);
    SDFgetint(sdfp2, "npart", &npart2);
    printf("%d and %d particles\n", npart1, npart2);


/*
	maxmbrs = nmembers1 > nmembers2 ? nmembers1 : nmembers2;
	minmbrs = nmembers1 > nmembers2 ? nmembers2 : nmembers1;
*/
	maxmbrs = nmembers1;
	minmbrs = nmembers2;
	printf("maxmbrs = %d\n", maxmbrs);


	/*check for correspondence between headers */
	sort_head = verify_head(vecs1, nvecs1, nmembers1, vecs2, nvecs2, nmembers2, &ordhead);

	/*check for correspondence between structs */
	sort_body = verify_body(vecs1, nvecs1, nmembers1, vecs2, nvecs2, nmembers2, &ordbody);

	if(nmembers1 == nmembers2)
		verify_abun(abarr1, abarr2, nmembers2);


    /* if nvecs1 != nvecs2 there are two options:
		- use shorter header and skip extra data
		- use longer header and fill extra data with zeros
	   For now, do option 2, as the extra data is bndry info that
	   I don't want to skip. -CIE
	*/


    /*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(maxmbrs * sizeof(char *));
    if(members == NULL) printf("allocation error\n");
    addrs = (void **)malloc(maxmbrs * sizeof(void *));
    types = (SDF_type_t *)malloc(minmbrs * sizeof(SDF_type_t));
    outtypes = (SDF_type_t *)malloc(maxmbrs * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(minmbrs * sizeof(int));
    outoffsets = (int *)malloc(maxmbrs * sizeof(int));
    strides = (int *)malloc(maxmbrs * sizeof(int));
    starts = (int *)malloc(maxmbrs * sizeof(int));
    lines = (int *)malloc(maxmbrs * sizeof(int));

    printf("done malloc'ing\n");

    flag = 0;
/* note-to-self: temporary solution, clean up eventually! ~CIE */
	/* prepare arguments for SDF-reading the appropriate file */
/*
    if ( (nmembers1 > nmembers2) || (nvecs1 > nvecs2) ) {
*/
		printf("doing file 1 first\n");
		second = 2;

		npart = npart1;
		sdfp = sdfp1;
		idindex = prep_SDF_args(sdfp, vecs1, nvecs1, &members, &nmembers,
								&outtypes, &outoffsets, &strides, &lines, &outstride);

/*
    } else {
		printf("doing file 2 first\n");
		second = 1;

		npart = npart2;
		sdfp = sdfp2;
		idindex = prep_SDF_args(sdfp, vecs2, nvecs2, &members, &nmembers,
								&outtypes, &outoffsets, &strides, &lines, &outstride);

    }
*/

    /* unnecesary, just use 'stride' ? CIE */
/*
    outstride = 0;
    for(i=0; i< maxmbrs; i++) outstride += SDFtype_sizes[ outtypes[i] ];
*/

    printf("outstride = %d\n", (int)outstride);

    /* holds the read in data */
    btab = (void *)malloc( outstride*incr );
    outbtab = (void *)malloc( outstride*incr );

    /*calculate the byte offset in memory to the address of the next member-CIE*/
    addrs[0] = (char *)btab;
    for (i=1; i< maxmbrs; i++)
         addrs[i] = addrs[i-1] + SDFtype_sizes[ outtypes[i-1] ];


    printf("printing header\n");
    writescalars(sdfp1, fp, &pos_npart,npart1+npart2);

    /* writes the new header from the appropriate file.
     * Returns the location in memory where the
     * total particle number is stored so it can be updated later */

/*
    if (nvecs2 > nvecs1) {
        writescalars(sdfp2, fp, &pos_npart,npart1+npart2);
    } else {
        writescalars(sdfp1, fp, &pos_npart,npart1+npart2);
    }
*/

    /*print the struct declaration part from the header-CIE*/
/* TO DO: update 'R0' */
    fprintf(fp, "struct {\n");
    for (i = 0; i < maxmbrs; ++i) {
	switch (outtypes[i]) {
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

    /*loop over first file in chunks of data, write that chunk to file*/
    printf("getting first file .... ");
    for( i=0, countnpart = 0; i < npart; i++) {
        /* need to increment starts-array */
        for( j = 0; j < maxmbrs; j++)
            starts[j] = i;

        /* read data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, maxmbrs, members, starts, lines, addrs, strides);

        ident = *((int *)(btab + outoffsets[idindex]));

        /* calculate radius for merging */
        x = *((double *)(btab + outoffsets[0]));
        y = *((double *)(btab + outoffsets[1]));
        z = *((double *)(btab + outoffsets[2]));
        radius = sqrt(x*x + y*y + z*z);

	/*dump the btab data into the file now-CIE*/
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

    printf("got %d lines\n",countnpart);

	SDFclose(sdfp1);


	/* prepare arguments for SDF-reading from second file.
	   Take special care that the 'vecs' are in the same order
	   as in the first file
	*/
	/* I just don't seem to have any luck with realloc ~CIE */
	free(starts);
	free(strides);
	free(lines);
	for( i = 0; i < maxmbrs; i++)
		free(members[i]);
	free(members);

	starts = (int *)malloc( minmbrs*sizeof(int) );
	strides = (int *)malloc( minmbrs*sizeof(int) );
	members = (char **)malloc( minmbrs*sizeof(char*) );
	lines = (int *)malloc( minmbrs*sizeof(int) );

/*
	if(second == 1) {
		printf("doing file 1 now\n");
		sdfp = sdfp1;
		npart = npart1;

		idindex = prep_SDF_args(sdfp, vecs1, nvecs1, &members, &nmembers,
								&types, &inoffsets, &strides, &lines, &stride);

		free(btab);
		btab = (void *)malloc(stride*incr );

*/
	    /*calculate the byte offset in memory to the address of the next member-CIE*/
/*
	    addrs[0] = (char *)btab;
	    for (i=1; i< minmbrs; i++)
	         addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];


	} else if( second == 2) {
*/
		printf("doing file 2 now\n");
		sdfp = sdfp2;
		npart = npart2;

		idindex = prep_SDF_args(sdfp, vecs2, nvecs2, &members, &nmembers,
								&types, &inoffsets, &strides, &lines, &stride);

		free(btab);
		btab = (void *)malloc(stride*incr );

	    /*calculate the byte offset in memory to the address of the next member-CIE*/
	    addrs[0] = (char *)btab;
	    for (i=1; i< minmbrs; i++)
	         addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];


/*
	} else printf("huh??\n");
*/

/* verify that nmembers == minmbrs*/


    printf("getting file 2 .... ");
    for( i=0, countnpart = 0; i < npart; i++) {
        /* need to increment starts-array */
        for( j = 0; j < nmembers; j++)
            starts[j] = i;

        /* read data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, maxmbrs, members, starts, lines, addrs, strides);

        /* update the particle id, so it does not start at 1 again */
        ident = *((int *)(btab + inoffsets[idindex]));
        ident += identmax + 1;
        memcpy( btab + inoffsets[idindex], &ident, sizeof(ident) );

        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        radius = sqrt(x*x + y*y + z*z);

		/*dump the btab data into the file now-CE*/
		/* sort to match first header first */
		if( sort_body ) {
			for( k=0; k < maxmbrs ; k++ ) {
				memcpy(outbtab + outoffsets[k], &sort_body, SDFtype_sizes[outtypes[k]]);
				if( ordbody[k] >= 0) {
					memcpy(outbtab + outoffsets[k], btab + inoffsets[ordbody[k]],
							SDFtype_sizes[types[ordbody[k]]]);
				} /*else {
					memcpy(outbtab + outoffsets[k], 0., types[ordbody]);
				}*/
			}
		} else {
			memcpy(outbtab, btab, outstride);
        }

        if( radius > R0 ) { /* radius is always greater than -99., no extra case needed*/
	    fwrite(outbtab, outstride, 1, fp);
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

	SDFclose(sdfp2);

/*and we're done! clean up now -CE: if it ever works*/
/*    free(members);
    free(btab);
    free(addrs);
*/
    free(types);
    free(inoffsets);
	free(starts);
	free(lines);
    /*free(outbtab);*/
	if(!ordhead) free(ordhead);
	if(!ordbody) free(ordbody);

}


int GetNumOfDigits(int number) {
    int digits = 1;

    while(number/10 >= 1) {
        number = number/10;
        digits++;
    }
    return digits;
}




int verify_head(char **vecs1, int nvecs1, int nmembers1, char **vecs2, int nvecs2, int nmembers2, int **ordhead)
{
	int i,j;


	if( (nvecs1-nmembers1) != (nvecs2-nmembers2) ) {

		printf("different content in headers! %d vs %d\n\n", nvecs1, nvecs2);
		for(i=0; i< (nvecs1 > nvecs2 ? nvecs2 : nvecs1); i++)
			printf("%d) %s   %s\n",i,vecs1[i], vecs2[i]);


		if( (nvecs1-nmembers1) > (nvecs2-nmembers2) ){

			j = 0;
			(*ordhead) = (int *)malloc( sizeof(int)*(nvecs1-nmembers1) );
			do{
				(*ordhead)[j] = -1;
				for( i=0; i< (nvecs2-nmembers2); i++) {
					if( strncmp(vecs1[j],vecs2[i],strlen(vecs1[j])) == 0)
						(*ordhead)[j]=i;
				}
				j++;
			} while(j < (nvecs1-nmembers1));

		} else {

			j = 0;
			(*ordhead) = (int *)malloc( sizeof(int)*(nvecs2-nmembers2) );
			do {
				(*ordhead)[j] = -1;
				for( i=0; i< (nvecs1-nmembers1); i++) {
					if( strncmp(vecs2[j],vecs1[i],strlen(vecs2[j])) == 0)
						(*ordhead)[j]=i;
				}
				j++;
			}while(j < (nvecs2-nmembers2));

		}

		return 1;

	} else {

	return 0;

	}

}


int verify_body(char **vecs1, int nvecs1, int nmembers1, char **vecs2, int nvecs2, int nmembers2, int **ordbody)
{

	/* at some point create a more rigorous check, i.e. for cases when the
	 * number of struct members is the same, but the content is not
	*/

	int off1, off2;
	int i, j, abundflag = 0;

	off1 = nvecs1 - nmembers1;
	off2 = nvecs2 - nmembers2;

	if( nmembers1 != nmembers2 ) {

		abundflag = 1;

		printf("different content in structs! %d vs %d\n\n", nvecs1, nvecs2);
		for(i=0; i< (nmembers1 > nmembers2 ? nmembers1 : nmembers2); i++)
			printf("%s   %s\n",vecs1[off1+i], vecs2[off2+i]);
/*
*/
	}

/*
		if( nmembers1 > nmembers2 ){
*/

			j = 0;
			(*ordbody) = (int *)malloc( sizeof(int)*nmembers1 );
			do {
				(*ordbody)[j] = -1;
				for( i=0; i< nmembers2; i++) {
					if( strncmp( vecs1[off1+j],vecs2[off2+i],
								strlen(vecs1[off1+j]) ) == 0)
						(*ordbody)[off2+i]=off1+j;
				}
				j++;
			}while(j < nmembers1) ;

/*
		} else {

			j = 0;
			(*ordbody) = (int *)malloc( sizeof(int)*nmembers2 );
			do {
				(*ordbody)[j] = -1;
				for( i=0; i< nmembers1; i++) {
					if( strncmp( vecs2[off2+j],vecs1[off1+i],
								strlen(vecs2[off2+j]) ) == 0)
						(*ordbody)[off1+i]=off2+j;
				}
				j++;
			}while(j < nmembers2) ;
*/

/*
		}
*/

		return 1;

/*
	} else {

		return 0;

	}
*/

}


int verify_abun(int abarr1[][22], int abarr2[][22], int nmembers)
{
	int i, k, abundflag = 0;


		for( i = 0; i < nmembers; i++ ) {
			if( (abarr1[0][i] != abarr2[0][i]) || (abarr1[1][i] != abarr2[1][i]) )
				abundflag = 1;
		}

        /* check if the abundance information is the same in both files, and spit
         * out a warning if it is not */
        if(abundflag == 1) {
           printf("warning: abundance info in files might not be the same!\n");
           printf("%10s %10s\n%5s%5s %5s%5s\n", "file1", "file2","p","n","p","n");
           for(k=0;k<22;k++) printf("%5d%5d %5d%5d\n",
               abarr1[0][k],abarr1[1][k],abarr2[0][k],abarr2[1][k]);
    }

	return 0;

}




int prep_SDF_args(SDF *sdfp, char **vecs, int nvecs, char ***members, int *nmembers, SDF_type_t **types, int **offsets, int **strides, int **lines, size_t *stride)
{

	int i, incr=1, flag = 0, idindex = 0;

        for (i = 0, *stride = 0, *nmembers = 0; i < nvecs; ++i) {

            if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag=1;

            if(flag) {
                (*members)[*nmembers] = vecs[i];
                (*types)[*nmembers] = SDFtype((*members)[*nmembers], sdfp);
                (*offsets)[*nmembers] = *stride;
                *stride += SDFtype_sizes[(*types)[*nmembers]];
                (*lines)[*nmembers] = incr;
                (*nmembers)++;
            }

            /* get index that holds 'ident' for later updating */
            if (strncmp(vecs[i], "ident", strlen(vecs[i])) == 0)
				idindex = *nmembers - 1;

        }
        for (i = 0; i < *nmembers; i++)
            (*strides)[i] = *stride;
		printf("stride = %d\n",(int)(*stride));

	return idindex;
}
