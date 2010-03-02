/*
   PURPOSE:
	do the below for a whole bunch of files.
	to convert SDF files to formatted text files that OpenDX can read.
        only a few columns are written to the text file, namely x,y,z,h,mass,rho.
	One line of the SDF file is read at a time, so this should work for any
	size SDF file.
   NOTE: this code assumes that 'npart' from the header contains the total number
	of particles.

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
#include <SDF.h>

typedef enum SDF_type_enum SDF_type_t;

typedef union {
    int i;
    float f;
    double d;
} datum_t;

double ***trajec;
int maxID = 200000, idstep=10000;

static void writestructs(SDF *sdfp, FILE *fp, int row);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;
    int i, j, k, start, last, INCR, *len;
    int nfiles, nrows;
    char ext[5],extID[8];
    char *sdfname,*outfile;


    /* check for correct arguments */
    if (argc != 3) {
	fprintf(stderr, "Usage: %s SDFfile outfile \n", argv[0]);
	exit(1);
    }

    len = (int *)malloc(argc*sizeof(int));
    len[1] = strlen(argv[1]);
    len[2] = strlen(argv[2]);
    printf("argv[1]: %s has length %d\n", argv[1], len[1]);

    //myid = 1000;
    INCR = 50;
    start = 0;
    last = 500;

    sdfname = (char *)malloc( (len[1]+5 ) * sizeof( char ) );
    outfile = (char *)malloc( (len[2]+8 ) * sizeof( char ) );

    strcpy(sdfname, "/scratch/cellinge/runsnsph/run3g_50a1_sph.");
    strcpy(outfile, "/scratch/cellinge/run3g_50.traj.");

    printf("%s  %d  %d\n", sdfname, strlen(sdfname),len[1]);
    printf("%s  %d  %d\n", outfile, strlen(outfile),len[2]);

    nfiles = (int)maxID/idstep;
    nrows = (int)(last-start)/INCR+1;
    printf("files: %d rows: %d\n", nfiles, nrows);

    trajec = (double ***)malloc( nfiles * sizeof( double **));
    for ( j = 0; j < nfiles; j++) {
         trajec[j] = (double **)malloc( nrows * sizeof( double * ) );
        for ( i = 0; i < nrows; i++)
             trajec[j][i] = (double *)malloc( 6 * sizeof(double) );
    }

    /* iterate over sdf files */
    for(i = 0; i <= last; i += INCR) {
        printf("%-4d: ",i);
        sprintf(&ext, "%04d", i);
        /*printf("ext: %s\n", ext);*/

        for( k = 0; k < 5; k++) {
            strncpy(&sdfname[ len[1]+k ], &ext[k], 1);
        }

        printf("opening: %s %d row: %d\n", sdfname,strlen(sdfname), (i/INCR));

        sdfp = SDFopen(NULL, sdfname);
        if (sdfp == NULL) {
	    fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
	    exit(2);
        }

        writestructs(sdfp, fp, (i/INCR));
        SDFclose(sdfp);
    }

    printf("all files read in\n");

    for( j = 0; j < nfiles; j++) {
        sprintf(&extID, "%-d", j*idstep);

        for( k = 0; k < 7; k++) strncpy(&outfile[ len[2]+k ], &extID[k], 1);
        k = 0;
        do {
           strncpy(&outfile[ len[2]+k ], &extID[k], 1);
           //printf("%s\n", extID[k]);
           k++;
        } while( (extID[k] != ' ') && (extID[k] != '\0') );
        //strncpy(&outfile[ len[2]+k ], '\0', 1);

        printf("\n\nwriting trajectory file %s\n", outfile);

        fp = fopen(outfile, "w");
        if(!fp) printf("error opening file!\n");
        fprintf(fp,"%13s %13s %13s %13s %13s %13s\n",
                "x",
                "y",
                "z",
                "rho",
                "temp",
                "time"
               );
        for( i = 0; i < nrows; i++) {
            printf("printing file %4d row %4d\n",j,i);
            fprintf(fp, "%+13E %+13E %+13E %+13E %+13E %+13E\n",
                    trajec[j][i][0],
                    trajec[j][i][1],
                    trajec[j][i][2],
                    trajec[j][i][3],
                    trajec[j][i][4],
                    trajec[j][i][5]
                   );
        }

    fclose(fp);
    }

    return 0;
}


/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp, int row)
{
    int i, j, nvecs, nmembers;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0, outstride = 0;
    void *btab;
    void **addrs;
    int *inoffsets, *lines, *strides, *starts;
    int INCR=1, flag=0, num=6;
    int nlines = 1, nrecs;
    int index[num];
    int partid;
    float timest, rho;
    int myid, fn;
    /*make INCR and nlines user input */
   /* printf("id: %d\n",myid);*/

/* this does not attempt to load the whole file into memory, does it? -CE */
    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    SDFgetfloat(sdfp, "tpos", &timest);
    timest = timest * 100.;
    printf("at time: %G\n", timest);

    /* Count structure members */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid. But I know that the structure (so far) always has "x" as the first
       member, so I can start counting from there -CE */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        /*get columns of interest*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) {
            /* x is the first member of the structure */
            index[0]=i;
            flag=1;
        }
        if (strncmp(vecs[i], "y", strlen(vecs[i])) == 0) index[1]=i;
        if (strncmp(vecs[i], "z", strlen(vecs[i])) == 0) index[2]=i;
        if (strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) index[3]=i;
        if (strncmp(vecs[i], "temp", strlen(vecs[i])) == 0) index[4]=i;
        if (strncmp(vecs[i], "ident", strlen(vecs[i])) == 0) index[5]=i;
	if (flag) ++nmembers;
    }

/*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(num * sizeof(char *));
    addrs = (void **)malloc(num * sizeof(void *));
    types = (SDF_type_t *)malloc(num * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(num * sizeof(int));
    strides = (int *)malloc(num * sizeof(int));
    lines = (int * )malloc(num * sizeof(int));
    starts = (int * )malloc(num * sizeof(int));

    /*printf("nmembers = %d\n",nmembers);*/
/*one by one, go through the fields in the column, i.e. members of the struct?-CE*/
    for (i = 0, stride = 0; i < num; ++i) {
	    members[i] = vecs[index[i]];
	    types[i] = SDFtype(members[i], sdfp);
	    inoffsets[i] = stride;/* offsets (from beginning of 'line'?) of each column
                                            of data (struct member) */
	    stride += SDFtype_sizes[types[i]];
            lines[i] = nlines;
	}

    /* unnecesary, just use 'stride' ? CE */
    outstride = 0;
    for( i = 0; i < num; i++) outstride += SDFtype_sizes[ types[i] ];

    btab = (void *)malloc(stride * nlines);

    /*calculate the byte offset in memory to the address of the next member-CE*/
	addrs[0] = (char *)btab;
	for ( i = 1; i < num; i++) {
            addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];
            strides[i] = outstride;
        }

    printf("reading in %d lines ...\n",nrecs);

    myid = 0;
    fn = 0;
    /*try reading in one line at a time */
    for( j = 0; j < nrecs; j++) {
        for( i = 0; i < num; i++)
            starts[i] = j;

	    SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);
            partid = *((int *)(btab + inoffsets[5]));
            //printf("%d  %E\n",partid,*(double *)(btab + inoffsets[0]));

            if(partid == myid) {

		trajec[fn][row][0] = *((double *)(btab + inoffsets[0]))*6.955e10;
		trajec[fn][row][1] = *((double *)(btab + inoffsets[1]))*6.955e10;
		trajec[fn][row][2] = *((double *)(btab + inoffsets[2]))*6.955e10;
		trajec[fn][row][3] = *((float *)(btab + inoffsets[3])) *
                       ((1.9889e27)/(6.955e10 * 6.955e10 * 6.955e10));
		trajec[fn][row][4] = *((float *)(btab + inoffsets[4]));
                trajec[fn][row][5] = timest;

                fn++;
                myid += idstep;
            }
            if(myid >= maxID) j = nrecs;
    }


/*and we're done! clean up now -CE: if it ever works*/
/*
    free(members);
    free(btab);
    free(addrs);
    free(types);
    free(inoffsets);
*/

}
