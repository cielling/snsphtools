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
int locate(float *bins, int Nbins, float val);


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

/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp)
{
    SDF_type_t *types, type;
    size_t stride = 0, outstride = 0;
    void *outbtab, *btab;
    void **addrs;
    int i, j, k, nvecs, nmembers;
    int INCR=1, flag=0, num=14;
    int Nbins = 100, bin;
    int *inoffsets, *lines, *strides, *starts;
    int nlines = 1, nrecs;
    int index[num];
    char **vecs, **members;
    char getmembrs[num][13];
    float **dm, *dv;
    float X_i[7], vbins[Nbins], mbins[Nbins][7], mtot[7];
    float vmin=100, vmax=0, delv;
    float vx, vy, vz;
    float vel,mass,masscf,lengthcf,timecf;
    double x, y, z;
    /*make INCR and nlines user input */

    for( k=0; k<(num-7); k++) mtot[k]=0.;

    /* specify here which quantities to read in.
       note: a few things downstairs depend on the
       order these are specified. change with care. */
    strcpy(getmembrs[0],"x");
    strcpy(getmembrs[1],"y");
    strcpy(getmembrs[2],"z");
    strcpy(getmembrs[3],"vx");
    strcpy(getmembrs[4],"vy");
    strcpy(getmembrs[5],"vz");
    strcpy(getmembrs[6],"mass");
    strcpy(getmembrs[7],"f19"); /* H */
    strcpy(getmembrs[8],"f2"); /* O */
    strcpy(getmembrs[9],"f3"); /* Si */
    strcpy(getmembrs[10],"f5"); /* Fe */
    strcpy(getmembrs[11],"f7"); /* Ni */
    strcpy(getmembrs[12],"f15"); /* Fe */
    strcpy(getmembrs[13],"f17"); /* Ni */


    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    SDFgetfloatOrDefault(sdfp, "massCF", &masscf, 1.9889e27);
    SDFgetfloatOrDefault(sdfp, "timeCF", &timecf, 1.e2);
    SDFgetfloatOrDefault(sdfp, "lengthCF", &lengthcf, 6.955e10);

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
        for( k = 1; k < num; k++) {
            if (strncmp(vecs[i], getmembrs[k], strlen(vecs[i])) == 0) {
                index[k]=i;
                break;
            }
        }
        if (flag) ++nmembers;

    }
    printf("nmembers = %d\n",nmembers);

/*malloc memory space for the respective features of the struct-CIE*/
    dm = (float **)malloc(nrecs*sizeof(float*));
    for( i = 0; i < nrecs; i++) dm[i] = (float *)malloc( (num-7) * sizeof( float ));

    dv = (float *)malloc(nrecs*sizeof(float));

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


    printf("reading in %d lines ...\n",nrecs);

    /*try reading in one line at a time */
    for( j = 0; j < nrecs; j++) {
        for( i = 0; i < num; i++)
            starts[i] = j;

        SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);

        /* calculate quantities to apply threshold, if any */
        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        vx = *((float *)(btab + inoffsets[3]));
        vy = *((float *)(btab + inoffsets[4]));
        vz = *((float *)(btab + inoffsets[5]));
        mass = *((float *)(btab + inoffsets[6]));
        X_i[0] = *((float *)(btab + inoffsets[7]));
        X_i[1] = *((float *)(btab + inoffsets[8]));
        X_i[2] = *((float *)(btab + inoffsets[9]));
        X_i[3] = *((float *)(btab + inoffsets[10]));
        X_i[4] = *((float *)(btab + inoffsets[11]));
        X_i[5] = *((float *)(btab + inoffsets[12]));
        X_i[6] = *((float *)(btab + inoffsets[13]));

        vel = (vx*x+vy*y+vz*z)/sqrt(x*x+y*y+z*z)*lengthcf/timecf*1.e-5;
        dv[j] = vel;

        if( vel < vmin) vmin = vel;
        if( vel > vmax) vmax = vel;

        /* running total of mass per element i */
        for( k = 0; k < (num-7); k++) {
            mtot[k] += mass*X_i[k];
            dm[j][k] = mass*X_i[k];
        }


    }

    delv = (vmax - vmin)/(float)(Nbins-1);
    printf("delta: %e, max: %e, min: %e\n", delv, vmax, vmin);

    for( k = 0; k<Nbins; k++) {
        vbins[k] = vmin+delv*(float)k;
        for( j = 0; j < (num-7); j++ )
            mbins[k][j] = 0.;
    }

    for( j = 0; j < nrecs; j++) {
        /* determine vel bin */
        bin = locate(vbins, Nbins, dv[j]);

        /* sum up the mass in vel bin */
        for( k = 0; k < (num-7); k++)
            mbins[bin][k] += dm[j][k]/(mtot[k]*dv[j]);
    }

    fprintf(fp, "%14s ", "vbins");
    for(k = 7; k < num; k++)
        fprintf(fp, "%14s ", getmembrs[k]);
    fprintf(fp, "\n");

    for( k = 0; k < Nbins; k++) {
        fprintf(fp, "%14e ", vbins[k]);
        for(j = 0; j < (num-7); j++)
            fprintf(fp, "%14e ", mbins[k][j]);
        fprintf(fp, "\n");
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




int locate(float *bins, int Nbins, float val) {
    int jl, jm, ju;

        /*quick bisection to locate the radial bins I'm in -CE */
        if (val >= bins[Nbins-1])
            jl = Nbins-1;
        else if (val <= bins[0])
            jl = 0;
        else {
            ju = Nbins - 1;
            jl = 0;
            while (ju-jl > 1) {
               jm = (ju+jl) >> 1;
               if ((val >= bins[jm]) == (bins[Nbins-1] >= bins[0]))
                  jl=jm;
               else
                  ju=jm;
            }
        }
    return jl;
}

