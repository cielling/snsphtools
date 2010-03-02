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

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp);
static void writestructs(SDF *sdfp, FILE *fp, int myid);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;
    int i, k, start, nfile, INCR, *len;
    int myid;
    char ext[5],input;
    char *sdfname;


    /* check for correct arguments */
    if (argc != 3) {
	fprintf(stderr, "Usage: %s SDFfile outfile \n", argv[0]);
	exit(1);
    }

    len = (int *)malloc(argc*sizeof(int));
    len[1] = strlen(argv[1]);
    printf("argv[1]: %s has length %d\n", argv[1], len[1]);

    /* open output file */
    if (access(argv[2], F_OK) == 0) {
        fprintf(stderr, "%s: %s exists; overwrite (y/n)? ", argv[0],
		argv[2]);
        input = getc(stdin);
        if ((input != 'y') && (input != 'Y'))
            exit(3);
    }

    fp = fopen(argv[2], "w");
    if (fp == NULL) {
	fprintf(stderr, "%s: %s\n", argv[2], strerror(errno));
	exit(errno);
    }

    myid = 1000;
    INCR = 50;
    start = 0;
    nfile = 5000;

    sdfname = (char *)malloc( (len[1]+5 ) * sizeof( char ) );
    strcpy(sdfname, "/scratch/cellinge/runsnsph/run3g_50a1_sph.");
    printf("%s  %d  %d\n", sdfname, strlen(sdfname),len[1]);

    /* iterate over sdf files */
    for(i = 0; i <= nfile; i += INCR) {
        sprintf(ext, "%04d", i);
        printf("ext: %s\n", ext);

        for( k = 0; k < 5; k++) {
            strncpy(&sdfname[ len[1]+k ], &ext[k], 1);
        }

        printf("opening: %s %d\n", sdfname,strlen(sdfname));

        sdfp = SDFopen(NULL, sdfname);
        if (sdfp == NULL) {
	    fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
	    exit(2);
        }
        writestructs(sdfp, fp, myid);
        SDFclose(sdfp);
    }

    fclose(fp);

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

/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp, int myid)
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
    float timest;
    static int first = 0;
    double x, y, z;
    float rho;
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

    printf("nmembers = %d\n",nmembers);
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

    if(!first) {
        fprintf(fp,"%13s %13s %13s %13s %13s %13s\n",
            members[0],
            members[1],
            members[2],
            members[3],
            members[4],
            "time"
           );
    }

    printf("reading in %d lines ...\n",nrecs);

    /*try reading in one line at a time */
    for( j = 0; j < nrecs; j++) {
        for( i = 0; i < num; i++)
            starts[i] = j;

	    SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);
            partid = *((int *)(btab + inoffsets[5]));
    //        printf("%d  %E\n",partid,*(double *)(btab + inoffsets[0]));

            if(partid == myid) {

		x = *((double *)(btab + inoffsets[0]))*6.955e10;
		y = *((double *)(btab + inoffsets[1]))*6.955e10;
		z = *((double *)(btab + inoffsets[2]))*6.955e10;
		rho = *((float *)(btab + inoffsets[3]))*
                       1.9889e27/(6.955e10 * 6.955e10 * 6.955e10);

		fprintf(fp,"%+13E %+13E %+13E %+13E %+13E %+13E\n",
                        x, y, z, rho,
			*((float *)(btab + inoffsets[4])),
			timest
                       );

                j = nrecs;
            }
    }

    first++;

/*and we're done! clean up now -CE: if it ever works*/
/*
    free(members);
    free(btab);
    free(addrs);
    free(types);
    free(inoffsets);
*/

}
