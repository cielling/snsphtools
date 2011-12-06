/*
   PURPOSE:
    this is the "batch" version of SDFtoASCII. You could probably write a 
    script that calls SDFtoASCII a bunch of times. Or you can do it this way....

   NOTE:
    - this code assumes that 'npart' from the header contains the total number
    of particles.
    - this routine can handle 4- or 5-digit file extensions to the sdf files, 
    but not both at the same time (yet?), so if you have more that 9999 time 
    steps, you have to call it twice, once for <9999, once for >9999.
    - does not have a switch to calculate log values of quantities
    - the program will prompt you for the first time step, the last time step
    (inclusive), and the time step increment to read in; e.g. entering '0 1000 10'
    at the prompt will read in the sdf dumps '..._sph.0000', '..._sph.0010', ...,
    '..._sph.1000'. 

   SYNTAX:
    SDFtoASCII-batch <sdfname-base> <outname-base>

   DONE: 1) get list of names in SDF file
   DONE: 2) read all scalars
   DONE: 3) read selected structure members, INCR lines at a time
   DONE: 4) calculate offsets, write adjusted members into buffer
   DONE: 5) loop over 3-4 until whole file is read (or seg fault is reached ;-P)
   DONE: 6) write scalars and buffer to file

   TO DO:
    - prompt user for quantities of interest to get (leave x, y, z, in)
    - make this independent of/scalable to arbitrary number of quantities
    - add log functionality
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
static void writestructs(SDF *sdfp, FILE *fp);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;
    int i, j, k, start, nfile, INCR, *len;
    char ext[5];
    char *outname, *sdfname;

    len = (int *)malloc(argc*sizeof(int));
    for( j = 1; j < argc; j++) {
        len[j-1] = strlen(argv[j]);
        printf("argv: %s length: %d\n", argv[j],len[j-1]);
    }
    //len[0] += 27;
    sdfname = (char *)malloc( (len[0]+5) * sizeof( char ) );
    outname = (char *)malloc( (len[1]+5) * sizeof( char ) );
    /*maybe I'm too stupid, but argv[j][i] didn't work -CE */
/* WHY IS THIS NOT WORKING???
    for( i = 0; i <= len[2]; i++) {
        printf("outname. i: %d\n", i);
        outname[i] = *(char *)(argv[2] + i);
        printf("%s %c\n", outname[i], outname[i]);
    }
    for( i = 0; i <= len[1]; i++) {
        printf("sdfname. i: %d\n", i);
        sdfname[i] = *(argv[1]++);
        printf("%s %c\n", sdfname[i], sdfname[i]);
    }
*/

/*
    strcpy(outname, "/scratch/cellinge/run3g_50Am6.dat.");
    strcpy(sdfname, "/scratch/cellinge/run3g_50Am6_sph.");
*/
    strcpy(outname, argv[2]);
    strcpy(sdfname, argv[1]);

    printf("enter numbers of: first-file last-file file-increment \n");
    scanf("%d %d %d", &start, &nfile, &INCR);
    printf("%d %d\n", INCR, nfile);

    for( i = start; i <= nfile; i = i+INCR ) {
        sprintf(ext, "%04d", i);
        printf("ext: %s\n", ext);

        for( k = 0; k < 5; k++) {
            strncpy(&sdfname[ len[0]+k ], &ext[k], 1);
            strncpy(&outname[ len[1]+k ], &ext[k], 1);
        }

        printf("strnpcy: %s %s\n", sdfname, outname);

        sdfp = SDFopen(NULL, sdfname);
        if (sdfp == NULL) {
	    fprintf(stderr, "%s: %s: %s\n", argv[0], sdfname, SDFerrstring);
	    exit(2);
        }

        fp = fopen(outname, "w");
        if (fp == NULL) {
	    fprintf(stderr, "%s: %s\n", outname, strerror(errno));
	    exit(errno);
        }
    /*initargs(argc, argv, &sdfp, &fp);*/

        writestructs(sdfp, fp);

        fclose(fp);
        SDFclose(sdfp);

    }

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
static void writestructs(SDF *sdfp, FILE *fp)
{
    int i, j, k, nvecs, nmembers;
    char **vecs, **members;
    SDF_type_t *types, type;
    size_t stride = 0, outstride = 0;
    void *outbtab, *btab;
    void **addrs;
    int *inoffsets, *lines, *strides, *starts;
    int INCR=1, flag=0, num=7;
    int nlines = 1, nrecs;
    int index[num];
    double rx, ry, rz;

    /* this does not attempt to load the whole file into memory, does it? -CE */
    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    nrecs=8192;

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
        if (strncmp(vecs[i], "vx", strlen(vecs[i])) == 0) index[3]=i;
        if (strncmp(vecs[i], "vy", strlen(vecs[i])) == 0) index[4]=i;
        if (strncmp(vecs[i], "vz", strlen(vecs[i])) == 0) index[5]=i;
        if (strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) index[6]=i;
	if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

    /*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(num * sizeof(char *));
    addrs = (void **)malloc(num * sizeof(void *));
    types = (SDF_type_t *)malloc(num * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(num * sizeof(int));
    strides = (int *)malloc(num * sizeof(int));
    lines = (int * )malloc(num * sizeof(int));
    starts = (int * )malloc(num * sizeof(int));

    /*one by one, go through the fields in the column, i.e. members of the struct?-CIE*/
    for (i = 0, stride = 0; i < num; ++i) {
	    members[i] = vecs[index[i]];
	    types[i] = SDFtype(members[i], sdfp);
	    inoffsets[i] = stride;/* offsets (from beginning of 'line'?) of each column
                                            of data (struct member) */
	    stride += SDFtype_sizes[types[i]];
            lines[i] = nlines;
	}

    /* unnecesary, just use 'stride' ? CIE */
    outstride = 0;
    for( i = 0; i < num; i++) outstride += SDFtype_sizes[ types[i] ];

    btab = (void *)malloc(stride * nlines);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
	addrs[0] = (char *)btab;
	for ( i = 1; i < num; i++) {
            addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];
            strides[i] = outstride;
        }

/*
    for ( k = 0; k < num; k++) fprintf(fp,"% 13s", members[k]);
    fprintf(fp, "\n");
*/

    printf("reading in %d lines ...\n",nrecs);

    /*try reading in one line at a time */
    for( j = 0; j < nrecs; j++) {
        for( i = 0; i < num; i++)
            starts[i] = j;

        SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);

/*
        rx = (*(double *)(btab + inoffsets[0]))*1.e3;
        ry = (*(double *)(btab + inoffsets[1]))*1.e3;
        rz = (*(double *)(btab + inoffsets[2]))*1.e3;
*/

        if( (rx >= 0.) && (ry >= 0.) && (rz >= 0.) ) {
        for( k = 0; k < num; k++) {
            type = SDFtype(members[k],sdfp);
            switch(type){
            case SDF_FLOAT:
                fprintf(fp," %+13E", *(float *)(btab + inoffsets[k]));
                break;
            case SDF_DOUBLE:
                fprintf(fp," %+13E", (*(double *)(btab + inoffsets[k])));
                break;
	        default:
	            fprintf(stderr, "%s: type not supported\n", type);
                exit(-1);
            }
        }
        fprintf(fp, "\n");
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
