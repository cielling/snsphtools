/*
   PURPOSE:
    to convert SDF files to formatted text files that OpenDX can read.
    only a few columns are written to the text file, namely x,y,z,h,mass,rho.
    One line of the SDF file is read at a time, so this should work for any
    size SDF file.

   NOTE:
    - this code assumes that 'npart' from the header contains the total number
    of particles.
    - the physical quantities retrieved from the sdf can be changed (however, 
    the first always has to be 'x'). If the number of variables retrieved is 
    changed, also change 'num' to reflect the new value. 
    - currently, the log-flag options are as follows
      * 0 = calculates linear abundances (i.e. just gets data from the SDF)
      * 1 = calculates the log of the abundances (meaning mass fractions)
      * 2 = converts the units to cgs for everything, but get mass fraction
      * 3 = calculates the mass density in cgs of the isotopes, leaves x,y,z
            in code units, though

   COMPILE:
    with Makefile2. un-comment appropriate lines.
    This routine needs some libraries from the tree code and SDF routines,
    so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
    once for serial use (i.e. without PAROS flag).

   RUN:
    SDFtoASCII <in-file.sdf> <out-file.sdf> log-flag

   DONE: 1) get list of names in SDF file
   DONE: 2) read all scalars
   DONE: 3) read selected structure members, INCR lines at a time
   DONE: 4) calculate offsets, write adjusted members into buffer
   DONE: 5) loop over 3-4 until whole file is read (or seg fault is reached ;-P)
   DONE: 6) write scalars and buffer to file
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
static void writestructs(SDF *sdfp, FILE *fp);

int logg;

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;

    initargs(argc, argv, &sdfp, &fp);

    writestructs(sdfp, fp);

    fclose(fp);
    SDFclose(sdfp);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp)
{
    char input;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s SDFfile outfile log_flag\n", argv[0]);
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

    logg = atoi(argv[3]);
    printf("%d ", logg);
    if(logg==0) printf("calculating linear abundances\n");
    if(logg==1) printf("calculating log abundances\n");
    if(logg==2) printf("calculating mass density in cgs\n");
    if(logg==3) printf("calculating mass density in cgs, leaving x,y,z in code units\n");

}

/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp)
{
    int i, j, k, nvecs, nmembers;
    char **vecs, **members;
    SDF_type_t *types, type;
    size_t stride = 0, outstride = 0;
    void *outbtab, *btab;
    void **addrs;
    int *inoffsets, *lines, *strides, *starts;
    int INCR=1, flag=0, num=11;
    float X_i[4], vbins[100], mbins[100], mtot[4];
    int Nbins = 100;
    float vmin, vmax;
    int nlines = 1, nrecs;
    int index[num];
    double x, y, z;
    float vel,mass,masscf,lengthcf,timecf;
    /*make INCR and nlines user input */

    for( k=0; k<4; k++) mtot[k]=0.;

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    /*
    SDFgetint(sdfp, "massCF", &masscf);
    SDFgetint(sdfp, "timeCF", &timecf);
    SDFgetint(sdfp, "lengthCF", &lengthcf);
    */
    masscf = 1.998e27;
    lengthcf = 6.955e10;
    timecf = 1.e2;

    printf("length= %e, mass= %e, time= %e\n",
        lengthcf, masscf, timecf);

    /* Count structure members */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid. But I know that the structure (so far) always has "x" as the first
       member, so I can start counting from there -CIE */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        /*get columns of interest*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) {
            /* x is the first member of the structure */
            index[0]=i;
            flag=1;
        }
                /* change the read in quantities here */
                /* if you change the number of quantities read in, change
                   'num' in the variable declarations also */
        if (strncmp(vecs[i], "y", strlen(vecs[i])) == 0) index[1]=i;
        if (strncmp(vecs[i], "z", strlen(vecs[i])) == 0) index[2]=i;
        if (strncmp(vecs[i], "vx", strlen(vecs[i])) == 0) index[3]=i;
        if (strncmp(vecs[i], "vy", strlen(vecs[i])) == 0) index[4]=i;
        if (strncmp(vecs[i], "vz", strlen(vecs[i])) == 0) index[5]=i;
        if (strncmp(vecs[i], "mass", strlen(vecs[i])) == 0) index[6]=i;
        if (strncmp(vecs[i], "f1", strlen(vecs[i])) == 0) index[7]=i; /* c */
        if (strncmp(vecs[i], "f2", strlen(vecs[i])) == 0) index[8]=i; /* O */
        if (strncmp(vecs[i], "f5", strlen(vecs[i])) == 0) index[9]=i; /* Si */
        if (strncmp(vecs[i], "f17", strlen(vecs[i])) == 0) index[10]=i; /* Ni */

        if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

/*malloc memory space for the respective features of the struct-CIE*/
    members = (char **)malloc(num * sizeof(char *));
    addrs = (void **)malloc(num * sizeof(void *));
    types = (SDF_type_t *)malloc(num * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(num * sizeof(int));
    strides = (int *)malloc(num * sizeof(int));
    lines = (int * )malloc(num * sizeof(int));
    starts = (int * )malloc(num * sizeof(int));

/*one by one, go through the fields in the column, i.e. members of the struct?-CIE*/
    for (i = 0, stride = 0; i < num; ++i) {
        members[i] = vecs[index[i]];        /* quantities of interest */
        types[i] = SDFtype(members[i], sdfp);        /* their data types */
        inoffsets[i] = stride;/* offsets (from beginning of 'line'?) of each column
                                    of data (struct member) */
        stride += SDFtype_sizes[types[i]];        /* strides in memory */
        lines[i] = nlines;
        printf("member = %s offset: %d\n",members[i],inoffsets[i]);
    }

    /* unnecesary, just use 'stride' ? CIE */
    outstride = 0;
    for( i = 0; i < num; i++) outstride += SDFtype_sizes[ types[i] ];

        /* make room for btab to hold the read in data */
    btab = (void *)malloc(stride * nlines);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
    addrs[0] = (char *)btab;
    strides[0] = outstride;                /* how did this run correctly without this line??? */
    for ( i = 1; i < num; i++) {
        addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];
        strides[i] = outstride;
    }

    for( k = 0; k < num; k++) {
        fprintf(fp,"%13s\t", members[k]);
    }
    fprintf(fp,"\n");

    printf("reading in %d lines ...\n",nrecs);

    /*try reading in one line at a time */
    for( j = 0; j < nrecs; j++) {
        for( i = 0; i < num; i++)
            starts[i] = j;

        SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);

        /* calculate quantities to apply threshold, if any */
        vx = *((double *)(btab + inoffsets[1]));
        vy = *((double *)(btab + inoffsets[2]));
        vz = *((double *)(btab + inoffsets[3]));
        mass = *((float *)(btab + inoffsets[4]));
        X_i[0] = *((float *)(btab + inoffsets[5]));
        X_i[1] = *((float *)(btab + inoffsets[6]));
        X_i[2] = *((float *)(btab + inoffsets[7]));
        X_i[3] = *((float *)(btab + inoffsets[8]));

        vel = (vx*x+vy*y+vz*z)/sqrt(x*x+y*y+z*z);        
        if( vel < vmin) vmin = vel;
        if( vel > vmax) vmax = vel;

        /* sum up the mass */
        for( k = 0; k < 4; k++) 
            mtot[k] += mass*X_i[k];

        /* determine vel bin */

        /* sum up the mass in vel bin */

        }

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
