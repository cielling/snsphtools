/*
   PURPOSE:
	to add an abundance table to an initial conditions sdf file.
        This version should also work for large files, since it reads in the
        sdf file line-by-line.

   DONE: 1) get list of columns in SDF file
   DONE: 2) read all scalars
   DONE: 3) read all structure members
   DONE: 4) read in abundances and radial bins
   DONE: 5) find correct radial bin
   DONE: 6) write scalars, sdf data and abundances to file
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

static const int NISO = 22;

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
    FILE *fap = NULL;
    int i, j, k, nvecs, nmembers, Amembers, ju, jm, jl;
    int numA, Nbins, len, counter, flag = 0;
    int nrecs;
    int *strides, *nobjs, *starts, *inoffsets, *outoffsets, *nparr, *nnarr;
    int want[2][NISO]={{0,1,2,6,7,8,10,12,14,15,16,18,20,20,22,24,26,26,26,27,28,28},/*Z*/
                       {1,0,2,6,7,8,10,12,14,16,16,18,20,24,22,24,26,30,32,29,28,30}}; /*A-Z*/
    char **vecs, **members, **outmembers;
    char mychar[5];
    SDF_type_t *types, *outtypes;
    size_t stride = 0, outstride = 0;
    void *btab, *outbtab;
    void **addrs;
    float **abundarr; /*should this be void * ? */
    double x, y, z, radius;
    float *radbin;

    SDFgetint(sdfp, "npart", &nrecs);
    printf("%d particles\n", nrecs);

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

    Amembers = 3*NISO;

/*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    outmembers = (char **)malloc(3 * NISO * sizeof(char *));
    for(j=0;j<3*NISO;j++) outmembers[j] = (char *)malloc(10 * sizeof(char));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    strides = (int *)malloc(nmembers * sizeof(int));
    nobjs = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    outtypes = (SDF_type_t *)malloc(3 * NISO * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));
    outoffsets = (int *)malloc(Amembers * sizeof(int));

/*one by one, go through all the fields in the column, i.e. members of the struct?-CE*/
    flag = 0;
    for (i = 0, nmembers = 0, stride = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) {
	    members[nmembers] = vecs[i];
	    nobjs[nmembers] = 1;/*nobjs[0] is the number of particles. are all elements
                                      in this nobjs array the same then??-CE:yes, but can be
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

    /*calculate the byte offset in memory to the address of the next member-CE*/
    for (i = 0; i < nmembers; ++i) {
	addrs[i] = (char *)btab + inoffsets[i];
	strides[i] = stride;
    }

/* up to here just the SDF file is read in, so everything stays the same -CE */

    /* create the extra columns for abundances 'n stuff */
    for (i=0; i < NISO; i++) {
        len = sprintf(mychar, "%s%-d", "f",(i+1));
        for(j=0;j<5;j++) outmembers[i][j]=mychar[j];
        outtypes[ i ] = SDF_FLOAT;

        len = sprintf(mychar, "%s%-d", "p",(i+1));
        for(j=0;j<5;j++) outmembers[i+NISO][j]=mychar[j];
        outtypes[ i+NISO ] = SDF_INT;

        len = sprintf(mychar, "%s%-d", "m",(i+1));
        for(j=0;j<5;j++) outmembers[i+NISO*2][j]=mychar[j];
        outtypes[ i+NISO*2 ] = SDF_INT;
    }

    outstride = stride;

    /*calculate at what byte-intervals the data should be written-CE*/
    /* adding enough columns for the abundances */
    for (i = 0; i < Amembers; ++i) {
	outoffsets[i] = outstride;
	outstride += SDFtype_sizes[ outtypes[ i ] ];
    }

    /*malloc enough space in memory for the array that holds the whole output data-CE*/
    printf("about to malloc %u bytes for outbtab\n", nrecs*outstride);
    outbtab = (void *)malloc( outstride );

    /* now need to read in the abundances, and other necessary stuff */

    /* get Z and N of isotopes - WORKS!!*/
    fap = fopen("abundancereadme", "r");
    if (!fap) printf("error opening file abundancereadme\n");

    /*read in first number, that's the number of isotopes in the file*/
    fscanf(fap, "%3d", &numA);

    nparr = (int *)malloc(numA * sizeof(int));
    nnarr = (int *)malloc(numA * sizeof(int));

    for ( i = 0; i < numA; i++) fscanf(fap, "%d", &nparr[i]);
    for ( i = 0; i < numA; i++) fscanf(fap, "%d", &nnarr[i]);

    fclose(fap);
    fap = NULL;


    /* read in the radial bins - WORKS!!*/
    fap = fopen("inputh.dat", "r");
    if (!fap) printf("error opening inputh.dat\n");

    /*malloc array and read in radial bins. don't know how many, but guess 1000
      and realloc later if we run out of space -CE */
    radbin = (float *)malloc(1000 * sizeof(float));

    Nbins = 0;
    do {
        fscanf(fap, "%21G", &radbin[Nbins++]);
        fscanf(fap, "%*21g"); /* reads in h's. not needed */
        if(!(Nbins % 1000)) realloc(radbin, (Nbins+1000)*sizeof(float));
    } while(radbin[Nbins-1] < 4.3e2);

    fclose(fap);
    fap = NULL;

    /* malloc 2D array for abundances[radial bin][isotope] -CE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    /* read in abundance data - WORKS!!*/
    fap = fopen("abun.dat", "r");
    if (!fap) printf("error opening file abun.dat\n");

    /*read in abundances, one set for each radial bin*/
    for (i=0; i< Nbins; i++){
        for (j=0; j < numA; j++) {
            fscanf(fap, "%12G", &abundarr[i][j]);
        }
    }

    fclose(fap);
    fap = NULL;

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
    for (i = 0; i < Amembers; ++i) {
	switch (outtypes[i]) {
	case SDF_INT:
	    fprintf(fp, "\tint %s;\n", outmembers[i]);
	    break;
	case SDF_FLOAT:
            fprintf(fp, "\tfloat %s;\n", outmembers[i]);
	    break;
	case SDF_DOUBLE:
	    fprintf(fp, "\tdouble %s;\n", outmembers[i]);
            break;
	default:
	    fprintf(stderr, "%s: type not supported\n", outmembers[i]);
	    exit(-1);
	}
    }
    fprintf(fp, "}[%d];\n", nobjs);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    for (j = 0; j < nrecs; ++j) {
        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
		     strides);

        /* calculate radius. this should work -CE */
        x = *(double *)(btab + inoffsets[0]);
        y = *(double *)(btab + inoffsets[1]);
        z = *(double *)(btab + inoffsets[2]);
        radius = sqrt( x*x + y*y + z*z );

        /*quick bisection to locate the radial bin I'm in -CE */
        if (radius >= radbin[Nbins-1])
            jl = Nbins;
        else if (radius <= radbin[0])
            jl = 0;
        else {
            ju = Nbins - 1;
            jl = 0;
            while (ju-jl > 1) {
               jm = (ju+jl) >> 1;
               if ((radius >= radbin[jm]) == (radbin[Nbins-1] >= radbin[0]))
                  jl=jm;
               else
                  ju=jm;
            }
        }

	for (i = 0; i < nmembers; ++i) {
            memcpy(outbtab + inoffsets[i], btab + inoffsets[i],
	           SDFtype_sizes[ types[i] ]);
            starts[i] = j+1;
	}
/* which index holds my bin number? -CE: jl */
        /* now populate outbtab with abundances 'n stuff */
        for (i = 0, counter = 0; i < numA; i++) {
            for (k = 0; k < NISO; k++) {
            if ( (nparr[i] == want[0][k]) &&
                 (nnarr[i] == want[1][k]) ) {

                /* fill in abundances */
                memcpy(outbtab + outoffsets[ k ],
                   &abundarr[jl][i] , SDFtype_sizes[ outtypes[ k ] ]);

                /* fill in nprotons in isotope */
                memcpy(outbtab + outoffsets[ k+NISO ],
                   &nparr[i] , SDFtype_sizes[ outtypes[ k+NISO ] ]);

                /* fill in nneutrons in isotope */
                memcpy(outbtab + outoffsets[ k+NISO*2 ],
                   &nnarr[i] , SDFtype_sizes[ outtypes[ k+NISO*2 ] ]);

                counter++;
            }
            }
        }
        /*dump the outbtab data into the file now-CE*/
        fwrite(outbtab, outstride, 1, fp);
    }
    printf("wrote %d abundances to outbtab\n", counter);

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
