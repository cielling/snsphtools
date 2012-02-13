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
float Lodders(int Z);
float get_solar(int Z, int N);

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
    if(logg==4) printf("calculating abundance relative to solar (Lodders)\n");

}

/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp)
{
    SDF_type_t *types, type;
    size_t stride = 0, outstride = 0;
    void *outbtab, *btab;
    void **addrs;
    int i, j, k, nvecs, nmembers;
    int *inoffsets, *lines, *strides, *starts;
    int INCR=1, flag=0, num=11;
    int nlines = 1, nrecs;
    int index[num], getZ[5], getN[5];
    char **vecs, **members;
    char getmembrs[num][12], npchr[5], nnchr[5];
    double x, y, z;
    float rho,mass,h,masscf,lengthcf,timecf;
    float sol_val, sn_val, sn_h;
    /*make INCR and nlines user input */


    /* specify here which quantities to read in.
       note: a few things downstairs depend on the
       order these are specified. change with care. */
    strcpy(getmembrs[0],"x");
    strcpy(getmembrs[1],"y");
    strcpy(getmembrs[2],"z");
    strcpy(getmembrs[3],"rho");
    strcpy(getmembrs[4],"mass");
    strcpy(getmembrs[5],"h");
    strcpy(getmembrs[6],"f2");
    strcpy(getmembrs[7],"f4");
    strcpy(getmembrs[8],"f5");
    strcpy(getmembrs[9],"f15");
    strcpy(getmembrs[10],"f19"); /* H */


    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    SDFgetfloatOrDefault(sdfp, "massCF", &masscf, 1.9889e27);
    SDFgetfloatOrDefault(sdfp, "timeCF", &timecf, 1.e2);
    SDFgetfloatOrDefault(sdfp, "lengthCF", &lengthcf, 6.955e10);

    printf("length= %e, mass= %e, time= %e\n",
        lengthcf, masscf, timecf);

            if(logg == 4) {
/* fancy loop to get the A of the isotope for calculating atomic mass of it */
                for(k = 6; k<num; k++) {
                    strcpy(npchr, getmembrs[k]);
                    npchr[0]='p';
                    strcpy(nnchr, getmembrs[k]);
                    nnchr[0]='n';
                    SDFgetint(sdfp, npchr, &getZ[k-6]);
                    SDFgetint(sdfp, nnchr, &getN[k-6]);
                    printf("isotope %s = %d, %d\n",getmembrs[k],getZ[k-6],getN[k-6]);
                }
            }
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
        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        rho = *((float *)(btab + inoffsets[3]));
        mass = *((float *)(btab + inoffsets[4]));
        h = *((float *)(btab + inoffsets[5]));

        /*convert to cgs*/
        if(logg==2 || logg==3) {
            x = x*lengthcf;
            y = y*lengthcf;
            z = z*lengthcf;
            rho = rho*masscf/(lengthcf*lengthcf*lengthcf);
/*
            mass = mass*masscf;
            h = h*lengthcf;
*/


            /*copy back to btab*/
            if(logg==2) {
                memcpy( btab+inoffsets[0], &x, sizeof(x));
                memcpy( btab+inoffsets[1], &y, sizeof(y));
                memcpy( btab+inoffsets[2], &z, sizeof(z));
                memcpy( btab+inoffsets[3], &rho, sizeof(rho));
/*
                memcpy( btab+inoffsets[4], &mass, sizeof(mass));
                memcpy( btab+inoffsets[5], &h, sizeof(h));
*/
            }
            if(logg==3)
                memcpy( btab+inoffsets[3], &rho, sizeof(rho));
        }
/*
                At this point different thresholds can be set for
                selecting particles, e.g.:
        if((x >= 0.0) && (y >= 0.0) && (z >= 0.0)) {
*/
        if( rho >= 0.0e-0) {
            for( k = 0; k < num; k++) {

                /* calculate log abundances if flag is set */
                if(logg==1 && (k>2)){
                    *((float *)(btab + inoffsets[k])) =
                        log10(*((float *)(btab + inoffsets[k]))+1.e-20);
                }

                if( (logg==2 || logg==3) && k>3) {
                    *((float *)(btab + inoffsets[k])) *= rho;
                }

                if( logg == 4 && k > 5 && k < num-1) {
                   sol_val = get_solar(getZ[k-6], getN[k-6]);
                   sn_val = *((float *)(btab + inoffsets[k]));
                   sn_h = *((float *)(btab + inoffsets[num-1]));
                   sn_val = sn_val/sn_h*(1./((float)(getZ[k-6]+getN[k-6])));
                   sn_val = sn_val/sol_val;
                   memcpy(btab + inoffsets[k], &sn_val, SDFtype_sizes[ types[k] ]);
                }

                type = SDFtype(members[k], sdfp);

                switch(type) {
                case SDF_FLOAT:
                    fprintf(fp, "%+13E\t", *(float *)(btab + inoffsets[k]));
                    break;
                case SDF_DOUBLE:
                    fprintf(fp, "%+13E\t", *(double *)(btab + inoffsets[k]));
                    break;
                default:
                    printf("no such type: %s\n", type);
                    exit(1);
                }
            }

            fprintf(fp,"\n");
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




float Lodders(int Z){

    int i, j;
    typedef struct{
        int Z, N;
        float astro, cosmo;
    } solar;
    solar data_Lod[100];

    for( i = 0; i < 100; i++ )
        data_Lod[i].Z = i+1;

    data_Lod[0].astro = 12.00;
    data_Lod[1].astro = 10.98;
    data_Lod[2].astro = 3.35;
    data_Lod[3].astro = 1.48;
    data_Lod[4].astro = 2.85;
    data_Lod[5].astro = 8.46;
    data_Lod[6].astro = 7.90;
    data_Lod[7].astro = 8.76;
    data_Lod[8].astro = 4.53;
    data_Lod[9].astro = 7.95;
    data_Lod[10].astro = 6.37;
    data_Lod[11].astro = 7.62;
    data_Lod[12].astro = 6.54;
    data_Lod[13].astro = 7.61;
    data_Lod[14].astro = 5.54;
    data_Lod[15].astro = 7.26;
    data_Lod[16].astro = 5.33;
    data_Lod[17].astro = 6.62;
    data_Lod[18].astro = 5.18;
    data_Lod[19].astro = 6.41;
    data_Lod[20].astro = 3.15;
    data_Lod[21].astro = 5.00;
    data_Lod[22].astro = 4.07;
    data_Lod[23].astro = 5.72;
    data_Lod[24].astro = 5.58;
    data_Lod[25].astro = 7.54;
    data_Lod[26].astro = 4.98;
    data_Lod[27].astro = 6.29;
    data_Lod[28].astro = 4.34;
    data_Lod[29].astro = 4.70;
    data_Lod[30].astro = 3.17;
    data_Lod[31].astro = 3.70;
    data_Lod[32].astro = 2.40;
    data_Lod[33].astro = 3.43;
    data_Lod[34].astro = 2.67;
    data_Lod[35].astro = 3.36;
    data_Lod[36].astro = 2.43;
    data_Lod[37].astro = 2.99;
    data_Lod[38].astro = 2.28;
    data_Lod[39].astro = 2.67;
    data_Lod[40].astro = 1.49;
    data_Lod[41].astro = 2.03;
    data_Lod[42].astro = 1.89;
    data_Lod[43].astro = 1.18;
    data_Lod[44].astro = 1.77;
    data_Lod[45].astro = 1.30;
    data_Lod[46].astro = 1.81;
    data_Lod[47].astro = 0.87;
    data_Lod[48].astro = 2.19;
    data_Lod[49].astro = 1.14;
    data_Lod[50].astro = 2.30;
    data_Lod[51].astro = 1.61;
    data_Lod[52].astro = 2.35;
    data_Lod[53].astro = 1.18;
    data_Lod[54].astro = 2.25;
    data_Lod[55].astro = 1.25;
    data_Lod[56].astro = 1.68;
    data_Lod[57].astro = 0.85;
    data_Lod[58].astro = 1.54;
    data_Lod[59].astro = 1.02;
    data_Lod[60].astro = 0.60;
    data_Lod[61].astro = 1.13;
    data_Lod[62].astro = 0.38;
    data_Lod[63].astro = 1.21;
    data_Lod[64].astro = 0.56;
    data_Lod[65].astro = 1.02;
    data_Lod[66].astro = 0.18;
    data_Lod[67].astro = 1.01;
    data_Lod[68].astro = 0.16;
    data_Lod[69].astro = 0.84;
    data_Lod[70].astro = -0.06;
    data_Lod[71].astro = 0.72;
    data_Lod[72].astro = 0.33;
    data_Lod[73].astro = 1.44;
    data_Lod[74].astro = 1.42;
    data_Lod[75].astro = 1.75;
    data_Lod[76].astro = 0.91;
    data_Lod[77].astro = 1.23;
    data_Lod[78].astro = 0.88;
    data_Lod[79].astro = 2.13;
    data_Lod[80].astro = 0.76;
    data_Lod[81].astro = 0.16;
    data_Lod[82].astro = -0.42;

    return data_Lod[Z-1].astro;
}



float get_solar(int Z, int N) {

    float solar_h, solar_a;

    solar_a = Lodders(Z);

    solar_h = pow(10., (solar_a-12.));

    return solar_h;

}

