/*

   PURPOSE:
	to convert from the old file layout of the sdf files (Z, N, abundance
	listed for each particle) to the new one (Z, N information is in the
	header). Moving the isotope information from the data part into the
	header of the sdf file decreases the file size by roughly half.

   COMPILE:
	with Makefile2. un-comment appropriate lines.
	This routine needs some libraries from the tree code and SDF routines,
	so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
	once for serial use (i.e. without PAROS flag).

   SYNTAX:
	convertff <infile.sdf> <outfile.sdf>

   METHOD:
	@ open <infile.sdf> and read in header
	@ read in isotope info from first line of data
	@ laboriously add the isotope info to the header
	@ write the rest of the data to <outfile.sdf>

   NOTE:
	assumes that npart from the header contains the actual number of particles.
	also, in the past the free'ing of allocated arrays has caused segfaults
    Should be working now (02-08-2011).
	Tested and was working for run3g_50Am6 (Jul-23-2011)

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


static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp);
static void writeinit(FILE *fp);
static void writescalars(SDF *sdfp, FILE *fp);
static void writestructs(SDF *sdfp, FILE *fp);
int make_spec_names(char ***specarr, char spec, int num);

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

/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp)
{
    int i, j, k, nvecs, nmembers, noutmembers;
    char **vecs, **members, **outmembers;
    SDF_type_t *types;
    size_t stride = 0, outstride = 0;
    void *btab, *outbtab;
    void **addrs;
    int *inoffsets, *outoffsets, *starts, *lines, *strides;
    int flag=0, INCR=1;
    int nlines, n_iso=0, *parr, *narr;
    char **pnames, **nnames;
	double x;

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nlines);

    /* Count structure members */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid. But I know that the structure (so far) always has "x" as the first
       member, so I can start counting from there -CIE */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        /*get columns of interest*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag=1;
        if ((strlen(vecs[i]) == 2) && (strncmp(vecs[i], "f", 1) == 0)) n_iso++;
        else if ((strlen(vecs[i]) == 3) && (strncmp(vecs[i], "f", 1) == 0)) n_iso++;
	if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);
    printf("n_iso = %d\n", n_iso);
    noutmembers = nmembers - 2*n_iso;
    printf("noutmembers = %d\n",noutmembers);

	/*malloc memory space for the respective features of the struct-CIE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    lines = (int *)malloc(nmembers * sizeof(int));
    strides = (int *)malloc(nmembers * sizeof(int));
    parr = (int *)malloc( n_iso * sizeof(int) );
    narr = (int *)malloc( n_iso * sizeof(int) );
    outmembers = (char **)malloc(noutmembers * sizeof(char *));
    outoffsets = (int *)malloc(noutmembers * sizeof(int));

	/*one by one, go through the fields in the column, i.e. members of the struct?-CIE*/
    flag = 0;
    for (i = 0, stride = 0, outstride = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) {
            /* x is the first member of the structure */
            flag=1;
        }
        if(flag) {
            members[nmembers] = vecs[i];
            lines[nmembers] = INCR;
            types[nmembers] = SDFtype(members[nmembers], sdfp);
            inoffsets[nmembers] = stride;/* offsets (from beginning of 'line'?) of each column
                                            of data (struct member) */
            stride += SDFtype_sizes[types[nmembers]];
            if(nmembers < noutmembers) {
               outmembers[nmembers] = vecs[i];
/*
               types[nmembers] = SDFtype(outmembers[nmembers], sdfp);
*/
               outoffsets[nmembers] = outstride;
               outstride += SDFtype_sizes[types[nmembers]];
            }
            ++nmembers;
        }
    }

    printf("outstride= %d\n", (int)outstride);

    printf("malloc'ing btab: %d\n", (int)stride);
    btab = (void *)malloc( stride*INCR );
    outbtab = (void *)malloc( outstride*INCR );

    /*calculate the byte offset in memory to the address of the next member-CIE*/
	for (i=0; i<nmembers; i++) {
	    addrs[i] = (char *)btab + inoffsets[i];
        strides[i] = stride;
    }

    printf("reading in %d lines ...\n",nlines);

    make_spec_names(&pnames, 'p', n_iso);
    make_spec_names(&nnames, 'n', n_iso);

    /*try reading in one line at a time */
    for(j = 0; j < nlines; j = j+INCR) {

        for( i = 0; i < nmembers; i++)
            starts[i] = j;

        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, lines, addrs,
		     strides);

    if( j == 0) {
         for( i = 0; i < n_iso; i++) {
             parr[i] = *(int *)(btab + inoffsets[ nmembers - 2*n_iso + i ]);
             narr[i] = *(int *)(btab + inoffsets[ nmembers - 1*n_iso + i ]);
             fprintf(fp, "int %s = %d;\n", pnames[i], parr[i]);
             fprintf(fp, "int %s = %d;\n", nnames[i], narr[i]);
             printf("%s  int %s = %d;\n", members[nmembers-2*n_iso+i],pnames[i], parr[i]);
             printf("%s  int %s = %d;\n", members[nmembers-1*n_iso+i],nnames[i], narr[i]);
        }

        /*print the struct declaration part from the header-CIE*/
        fprintf(fp, "struct {\n");
        for (i = 0; i < noutmembers; ++i) {
            switch (types[i]) {
                case SDF_INT:
                    fprintf(fp, "\tint %s;\n", outmembers[i]);
	            printf("\tint %s;\n", outmembers[i]);
	            break;
	        case SDF_FLOAT:
                    fprintf(fp, "\tfloat %s;\n", outmembers[i]);
                    printf("\tfloat %s;\n", outmembers[i]);
	            break;
	        case SDF_DOUBLE:
	            fprintf(fp, "\tdouble %s;\n", outmembers[i]);
	            printf("\tdouble %s;\n", outmembers[i]);
	            break;
	        default:
	            fprintf(stderr, "%s: type not supported for member %d\n", outmembers[i], i);
	            exit(-1);
            }
	}
        /*print the struct declaration part from the header-CE*/
        fprintf(fp, "}[%d];\n", nlines);
        fprintf(fp, "#\n");
        fprintf(fp, "# SDF-EOH\n");

	printf("made p/n specifier\n");
        }

        for( k = 0; k < noutmembers; k++) {
            memcpy(outbtab + outoffsets[k], btab+inoffsets[k],
               SDFtype_sizes[ types[k] ]);
        }

        /*dump the outbtab data into the file now-CE*/
        fwrite(outbtab, outstride, INCR, fp);

    }


/*and we're done! clean up now -CE: if it ever works*/
/*
    for (i = 0; i < nmembers; i++) {
		free(members[i]);
		free(addrs[i]);
	}
    for (i = 0; i < noutmembers; i++) {
		free(outmembers[i]);
        }
*/
    printf("something doesnt wanna be free'd ....\n");
    free(members);
    free(addrs);
    free(types);
    free(inoffsets);
    free(starts);
    free(lines);
    free(strides);
	free(parr);
	free(narr);
	free(outoffsets);
    free(btab);
	free(outbtab);

}


int make_spec_names(char ***specarr, char spec, int num)
{
    int i;
    char tmpchr[20];

    *specarr = (char **)malloc(num * sizeof(char *) );

    for( i = 0; i < num; i++ ){

        sprintf( tmpchr, "%c%-d\0", spec, (i+1));
        *(*specarr + i) = (char *)malloc( ( strlen(tmpchr) + 1) * sizeof(char) );
        sprintf((*specarr)[i], "%s",tmpchr);
        //printf("isotope specifier: %s  %d\n", specarr[i],(int)strlen(specarr[i]));

    }
    return i;
}
