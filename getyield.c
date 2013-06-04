/*
   PURPOSE:
    to get the yield of a list of isotopes, and to calculate total, post-processed
    yields, if present

   NOTE:
    - this code assumes that 'npart' from the header contains the total number
    of particles.
    - the <isotopes-file> is just a table with the Z and N of each isotope to
    get the yields for
    - the <post-processed-file> is a file with the post-processed yields as
    currently given by the se_yields code. This is essentially a header that
    contains information about which particles were post-processed, and the yield
    for each isotope that was post-processed (I think). To get the un-Burn-ed
    yields, create a dummy post-processed file with two lines; the first line
    containing some (non-whitespace) characters, the second line containing just
    'nn'.

   COMPILE:
    make ARCH=<arch-value> [CC=<c-compiler>] PROGS=getyield -f Makefile

    - make sure $TREEHOME is set correctly, and that the tree library in that
    path has been compiled for serial use.

   SYNTAX:
    getyield <sdf-file> <output-file> <isotopes-file> <post-processed-file>

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

static int **isotope;
static int countis;

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp, FILE **fp2);
static void writestructs(SDF *sdfp, FILE *fp, FILE *fp2, int *partids, int npartids, float *ppyields);

int read_se_yields(char *argv[], int **partids, float **ppyields, FILE *fp2);
int make_spec_names(char ***specarr, char spec, int num);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL, *fp2 = NULL;
    int *partids, npartids;
    float *ppyields;

    initargs(argc, argv, &sdfp, &fp, &fp2);

    npartids = read_se_yields(argv, &partids, &ppyields, fp2);
    if( npartids <= 0)
        printf("WARNING! no post-proc'd particles found!\n");

    fclose(fp2);
    fp2 = fopen(argv[3], "r");

    writestructs(sdfp, fp, fp2, partids, npartids, ppyields);

    fclose(fp);
    fclose(fp2);
    SDFclose(sdfp);

    free(partids);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp, FILE **fp2)
{
    char input;

    if (argc != 5) {
        fprintf(stderr, "Usage: %s SDFfile outfile ctl-file postproc-file\n", argv[0]);
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

    *fp2 = fopen(argv[3], "r");
    if (*fp2 == NULL) {
        fprintf(stderr, "%s: %s\n", argv[3], strerror(errno));
        exit(errno);
    }

}

int read_se_yields(char *argv[], int **partid, float **ppyields, FILE *fp2)
{
    FILE *fp=NULL;
    int i, j, nlines=0;
    char temp[21], line[101], *lastpart;
    int counti=0, nn, np;
    float mass;


    /* read get.isotopes-ctl file */
    while( !feof(fp2) ){
        fgets(line, 100, fp2);
        counti++;
    }
    rewind(fp2);
    counti--;

    isotope = (int **)malloc( 2*sizeof(int *) );
    isotope[0] = (int *)malloc( counti*sizeof(int) );
    isotope[1] = (int *)malloc( counti*sizeof(int) );

    *ppyields = (float *)calloc( counti,sizeof(float) );


    counti = 0;

    while( !feof(fp2) ){
        fscanf(fp2, "%d %d", &isotope[0][counti], &isotope[1][counti]);
        counti++;
    }
    counti--;

    rewind(fp2);


    /* read post-processed particles file */
    fp = fopen(argv[4], "r");
    if(fp == NULL) {
        printf("error opening post-processing file! %s\n",argv[4]);
        return 0;
    }


    /* get IDs of particles that were post-processed */
    while( !feof(fp) ) {
        fgets(line, 100, fp);
        nlines++;
        lastpart = strstr(line, "nn");
        if(lastpart != NULL) break;
    }

    rewind(fp);
    nlines -= 2; /*step one back from '++', and one from header line */

    *partid = (int *)malloc( nlines*sizeof(int) );
    fgets(line, 100, fp); /* read in header line */

    for(i = 0; i < nlines; i++)
        fscanf(fp, "%s %d %*s %*s %*e", temp, &((*partid)[i]));


    /* get (total) post-processed masses */
    while( !feof(fp) ) {
        fscanf(fp, "%*s%*s%d%*s%*s%d%*s%*s%e%*s", &nn, &np, &mass);
        for( j = 0; j < counti; j++) {
            if((nn == isotope[1][j]) && (np == isotope[0][j])) {
            (*ppyields)[j] = mass/1.9889e33;
/*
                printf("mass= %e, %d %d\n",(*ppyields)[j], np, nn);
*/
            }
        }
    }


    fclose(fp);

    return nlines;
}


/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp, FILE *fp2, int *partids, int npartids, float *ppyields)
{
    int i, j, k, nvecs, nmembers;
    int ii, postprocd;
    char **vecs, **members;
    char **pnames, **nnames;
    char line[101];
    SDF_type_t *types, type;
    size_t instride = 0, outstride = 0;
    void *btab;
    void **addrs;
    int *inoffsets, *lines, *strides, *starts;
    int INCR=1, flag=0;
    int nlines = 1, nrecs;
    int index[6], counti, countlines, ident, ndx;
    int n_iso, **SDF_iso_arr, *iso_index;
    float rho, *yield, mass, *Xel, massCF, temp, total;
    float *oldyield;
    int **isotop;
    /*make INCR and nlines user input */

/* this does not attempt to load the whole file into memory, does it? -CE */
    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        /*get columns of interest*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) {
            /* x is the first member of the structure */
            index[0]=i;
            flag=1;
        }
        if (strncmp(vecs[i], "mass", strlen(vecs[i])) == 0) index[1]=i;
        if (strncmp(vecs[i], "ident", strlen(vecs[i])) == 0) index[2]=i;
        if (strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) index[3]=i;
        if (strncmp(vecs[i], "f1", strlen(vecs[i])) == 0) index[4]=i;
        if (strncmp(vecs[i], "p1", strlen(vecs[i])) == 0) index[5]=i;

        if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

    /* check for old file format vs new file format. In the old file
       format, Z and N of the isotopes was in the body of the data (i.e.
       members of the SPHbody struct, and followed the mass fraction f.
       In the new file format, they were added to the header, just before
       the struct.
    */
    if(index[5] > index[0]) { /* old format, p,n in body of struct */
        n_iso = (index[5]-index[4])/2;
        make_spec_names(&nnames, 'm', n_iso);
        make_spec_names(&pnames, 'p', n_iso);
    } else { /* new format */
        n_iso = (index[0]-index[5])/2;
        make_spec_names(&nnames, 'n', n_iso);
        make_spec_names(&pnames, 'p', n_iso);
    }

    //n_iso = 20;
    fprintf(stdout, "%d\n", n_iso);
    /* make space to hold the Z and N of each isotope from the SDF */
    SDF_iso_arr = (int**)malloc( 2*sizeof(int *) );
    SDF_iso_arr[0] = (int *)malloc( n_iso*sizeof(int) );
    SDF_iso_arr[1] = (int *)malloc( n_iso*sizeof(int) );


    if( SDFgetfloat(sdfp, "massCF", &massCF) != 0) massCF = 1.9889e27;
    massCF = massCF/1.9889e33;
    printf("massCF= %e\n", massCF);

    for( i=0; i < n_iso; i++) {
        SDFgetint(sdfp,pnames[i], &SDF_iso_arr[0][i]);
        SDFgetint(sdfp,nnames[i], &SDF_iso_arr[1][i]);
        fprintf(stdout, "%3s = %2d\t", pnames[i], SDF_iso_arr[0][i]);
        fprintf(stdout, "%3s = %2d\n", nnames[i], SDF_iso_arr[1][i]);
    /* at some point, it would be cool to print out the chemical symbol
       that corresponds to the given (np,nn) ... */
    }

    countlines = 0;
    while( !feof(fp2) ){
        fgets(line, 100, fp2);
        countlines++;
    }
    rewind(fp2);

    printf("getting %d isotopes\n",--countlines);
    //counti = countis+1;

    isotop = (int **)malloc( 2*sizeof(int *) );
    isotop[0] = (int *)malloc( countlines*sizeof(int) );
    isotop[1] = (int *)malloc( countlines*sizeof(int) );

    iso_index = (int *)malloc((countlines+3)*sizeof(int));
    Xel = (float *)malloc( countlines*sizeof(float) );
    yield = (float *)malloc( countlines*sizeof(float) );
    oldyield = (float *)malloc( countlines*sizeof(float) );

    for( i = 0; i < 3; ++i )
        iso_index[i] = index[i+1];

    counti=3; /* first three are mass, ident, rho indices */
    ii = 0;

    while( !feof(fp2) ){
        fscanf(fp2, "%d %d", &isotop[0][ii], &isotop[1][ii]);
        //fscanf(fp2, "%*d %*d");

        for(k = 0; k < n_iso; k++) {
            /* find isotope index in the SDF struct */
            if( (SDF_iso_arr[0][k] == isotop[0][ii]) &&
                (SDF_iso_arr[1][k] == isotop[1][ii]) ) {
                iso_index[counti++] = k + index[4];
                break;
            }
        }
        ii++;
    }

    fprintf(stdout, "counti: %d\n",counti);

    /* bubble- sort, to put isotope indices in same order as in the SDF */
    if(counti>3) {
        for(i = 3; i < counti; i++) {
            for(j = 3; j < counti-1; j++ ) {
                if(iso_index[j] > iso_index[j+1]) {
                    temp = iso_index[j+1];
                    iso_index[j+1] = iso_index[j];
                    iso_index[j] = temp;
                }
            }
        }
    }


    /*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc((counti) * sizeof(char *));
    addrs = (void **)malloc((counti)* sizeof(void *));
    types = (SDF_type_t *)malloc((counti) * sizeof(SDF_type_t));
    inoffsets = (int *)malloc((counti)* sizeof(int));
    strides = (int *)malloc((counti) * sizeof(int));
    lines = (int * )malloc((counti) * sizeof(int));
    starts = (int * )malloc((counti) * sizeof(int));

/*one by one, go through the fields in the column, i.e. members of the struct?-CE*/

    flag = 0;
    for (nmembers = 0, outstride = 0, j = 0; j < counti; ++j) {
        //members[j] = (char*)malloc( (strlen(vecs[iso_index[0]])+1)*sizeof(char) );
        //strcpy( members[j], vecs[iso_index[j]] );
        members[j] = vecs[iso_index[j]]; /* appears to be a memory leak. why? */
        lines[j] = 1;//INCR;
        starts[j] = 0;
        types[j] = SDFtype(members[j], sdfp);
        inoffsets[j] = outstride;
        outstride += SDFtype_sizes[ SDFtype(members[j],sdfp) ];
        fprintf(stdout, "member = %s offset: %d \n",members[j],inoffsets[j]);
        nmembers++;
    }

    fprintf(stdout, "nmembers: %d\n", nmembers);
    /* get the stride in bytes between successive entries in the same column */
    //for( i=0, instride = 0; i<nvecs; ++i )
    //    instride += SDFtype_sizes[ SDFtype(vecs[i], sdfp) ];
    instride = outstride;

    for( i=0; i<countlines; ++i) {
        /* initial yield-array to zero, since it carries a cumulative total */
        yield[i]=0.;
        oldyield[i] = 0.;
    }

    btab = (void *)malloc(outstride * nlines);

    /*calculate the byte offset in memory to the address of the next member-CE*/
    addrs[0] = (char *)btab;
    strides[0] = instride;
    for ( i = 1; i < nmembers; i++) {
        addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];
        strides[i] = instride;
    }


    printf("reading in %d lines ...\n",nrecs);


    total = 0.;

    /*try reading in one line at a time */
    for(ii = 0, j = 0; j < nrecs; j++) {
        for( i = 0; i < nmembers; i++)
            starts[i] = j;

        if(!(j%50000)) printf(".");
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, lines, addrs, strides);
            /*temporary solution*/

        mass = *((float *)(btab + inoffsets[0]));
        ident = *((int *)(btab + inoffsets[1]));
        rho = *((float *)(btab + inoffsets[2]));

        total += mass;

        postprocd = 0; /* only add abundances of non-post proc'd particles */

        if( npartids > 0 ) {
            if( (partids[ii] == ident) && (ii < npartids)) {
                postprocd = 1; /* don't add */
                ii++;
            } else {
                postprocd = 0; /* add */
            }
        }

        if( postprocd==0 ) {
            for(k = 3; k< counti; k++) {
                Xel[k-3] = *((float *)(btab + inoffsets[k]) );
                yield[k-3] += Xel[k-3] * mass * massCF;
            }
        } else {
            for( k=3; k<counti; k++ ) {
                Xel[k-3] = *((float *)(btab + inoffsets[k]) );
                oldyield[k-3] += Xel[k-3] * mass * massCF;
            }
        }

    }

    fprintf(fp, "%4s %4s %14s\t%14s\t%14s\t%14s\n", "p", "n", "un-proc'd", "proc'd", "orig", "total");
        for(k = 0; k < counti-3; k++) {
            ndx = iso_index[k+3]-index[4]; /* b/c we sorted this upstairs */
            fprintf(fp, "%4d %4d %14.6e", SDF_iso_arr[0][ ndx ],
                SDF_iso_arr[1][ ndx ], yield[k]);
            for(ii = 0; ii < countlines; ii++) {
                if((SDF_iso_arr[0][ ndx ] == isotop[0][ii]) &&
                    (SDF_iso_arr[1][ ndx ] == isotop[1][ii]) ) {
                    fprintf(fp, "\t%14.6e\t%14.6e\t%14.6e\n", ppyields[ii], oldyield[k], (yield[k]+ppyields[ii]));
                }
            }
        }
        fprintf(fp, "\n\ntotal mass: %e\n",total);

        printf("total mass: %e\n",total);

/*and we're done! clean up now -CE: if it ever works*/
/*    free(members);
    free(btab);
    free(addrs);
*/
    /*free(outbtab);*/
    free(lines);
    free(types);
    free(inoffsets);
    free(strides);
    free(SDF_iso_arr[0]);
    free(SDF_iso_arr[1]);
    free(SDF_iso_arr);
    //free(yield);
    //free(Xel);
    //free(oldyield);
    //free(pnames);
    //free(nnames);
    //free(ppyields);
    //free(isotope[1]);
    //free(isotope[0]);
    //free(isotope);

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
