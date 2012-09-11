/*
   PURPOSE:
    to create data files for the Digistar Planetarium software to render 
    visualizations/ movies of my simulations. Right now, the file format is
    the .vla format

   NOTE:
    - this code assumes that 'npart' from the header contains the total number
    of particles.
    - this routine can handle 4- or 5-digit file extensions to the sdf files, 
    but not both at the same time (yet?), so if you have more that 9999 time 
    steps, you have to call it twice, once for <9999, once for >9999.
    - the program will prompt you for the first time step, the last time step
    (inclusive), and the time step increment to read in; e.g. entering '0 1000 10'
    at the prompt will read in the sdf dumps '..._sph.0000', '..._sph.0010', ...,
    '..._sph.1000'. 

   SYNTAX:
    SDFtoASCII-batch <sdfname-base> <outname-base>

   DONE: 1) get list of names in SDF file
   DONE: 2) read all scalars
   DONE: 3) read selected structure members, INCR lines at a time
         3a) compute the "delta vector": r_[i+i] - r_[i]
   DONE: 4) calculate offsets, write adjusted members into buffer
   DONE: 5) loop over 3-4 until whole file is read (or seg fault is reached ;-P)
   DONE: 6) write scalars and buffer to file
         6a) write .vla file as 
               current x,y,z
               next x,y,z (or current x,y,z + delta x,y,z)
         6b) since Digistar's coord. sys. is left-handed, need to flip the sign 
             on y-components

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

static int first_file;

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp);
static void writestructs(SDF *sdfp, SDF *sdfp1, FILE *fp);
static void writeheader(FILE *fp);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL, *sdfp1 = NULL;
    FILE *fp = NULL;
    int i, j, k, start, nfile, INCR, *len;
    char ext[5];
    char *outname, *sdfname;

    sdfname = (char *)malloc( (strlen(argv[1])+10) * sizeof( char ) );
    outname = (char *)malloc( (strlen(argv[2])+10) * sizeof( char ) );
/*
    strcpy(outname, argv[2]);
    strcpy(sdfname, argv[1]);
*/

    printf("enter numbers of: first-file last-file file-increment \n");
    scanf("%d %d %d", &start, &nfile, &INCR);
    printf("%d %d\n", INCR, nfile);

    for( i = start; i <= nfile; i = i+INCR ) {
        if(i == start)
            first_file = 1;
        else
            first_file = 0;

        sprintf(sdfname, "%s%04d", argv[1], i);
        sprintf(outname, "%s%04d.vlc", argv[2], i);


        /* current position */
        sdfp = SDFopen(NULL, sdfname);
        if (sdfp == NULL) {
	    fprintf(stderr, "%s: %s: %s\n", argv[0], sdfname, SDFerrstring);
	    exit(2);
        }


        sprintf(sdfname, "%s%04d", argv[1], i+INCR);

        printf("strnpcy: %s %s\n", sdfname, outname);

        /* next position */
        sdfp1 = SDFopen(NULL, sdfname);
        if (sdfp1 == NULL) {
	    fprintf(stderr, "%s: %s: %s\n", argv[0], sdfname, SDFerrstring);
	    exit(2);
        }

        fp = fopen(outname, "w");
        if (fp == NULL) {
	    fprintf(stderr, "%s: %s\n", outname, strerror(errno));
	    exit(errno);
        }
    /*initargs(argc, argv, &sdfp, &fp);*/

        writeheader(fp);

        writestructs(sdfp, sdfp1, fp);

        fclose(fp);
        SDFclose(sdfp);
        SDFclose(sdfp1);

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
static void writestructs(SDF *sdfp, SDF *sdfp1, FILE *fp)
{
    int i, j, k, l, nvecs, nmembers;
    char **vecs, **members;
    SDF_type_t *types, type;
    size_t stride = 0, outstride = 0, instride = 0;
    void *outbtab, *btab, *btab1;
    void **addrs, **addrs1;
    int *inoffsets, *lines, *strides, *starts, *starts1;
    int INCR=1, flag=0, num=8;
    int nlines = 1, nrecs, nrecs1;
    int ident, ident1;
    int index[num];
    static int *idents = NULL;
    double x, x1, y, y1, z, z1, dx, dy, dz;
    float vx, vy, vz, vx1, vy1, vz1, dvel;
    float offs_x, offs_y, offs_z;

    printf("first_file is %d\n", first_file);

    /* this does not attempt to load the whole file into memory, does it? -CE */
    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    SDFgetint(sdfp1, "npart", &nrecs1);

    if( !idents)
        idents = (int *)calloc(nrecs,sizeof(int));

    /* Count structure members */
    /* don't use SDFnrecs, since that reads in the entire file which I'm trying to
       avoid. But I know that the structure (so far) always has "x" as the first
       member, so I can start counting from there -CE */
    for (i = 0, nmembers = 0, instride = 0; i < nvecs; ++i) {
        /*get columns of interest*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) {
            /* x is the first member of the structure */
            index[0]=i;
            flag=1;
        }
	if (flag) {
            ++nmembers;
            instride += SDFtype_sizes[SDFtype(vecs[i], sdfp)];
        }
        if (strncmp(vecs[i], "y", strlen(vecs[i])) == 0) index[1]=i;
        if (strncmp(vecs[i], "z", strlen(vecs[i])) == 0) index[2]=i;
        if (strncmp(vecs[i], "vx", strlen(vecs[i])) == 0) index[3]=i;
        if (strncmp(vecs[i], "vy", strlen(vecs[i])) == 0) index[4]=i;
        if (strncmp(vecs[i], "vz", strlen(vecs[i])) == 0) index[5]=i;
        if (strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) index[6]=i;
        if (strncmp(vecs[i], "ident", strlen(vecs[i])) == 0) index[7]=i;
    }
    printf("nmembers = %d\n",nmembers);

    /*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(num * sizeof(char *));
    addrs = (void **)malloc(num * sizeof(void *));
    addrs1 = (void **)malloc(num * sizeof(void *));
    types = (SDF_type_t *)malloc(num * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(num * sizeof(int));
    strides = (int *)malloc(num * sizeof(int));
    lines = (int * )malloc(num * sizeof(int));
    starts = (int * )malloc(num * sizeof(int));
    starts1 = (int * )malloc(num * sizeof(int));

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
    btab1 = (void *)malloc(stride * nlines);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
	addrs[0] = (char *)btab;
        strides[0] = outstride;
	for ( i = 1; i < num; i++) {
            addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];
            strides[i] = outstride;
        }
	addrs1[0] = (char *)btab1;
	for ( i = 1; i < num; i++) {
            addrs1[i] = addrs1[i-1] + SDFtype_sizes[ types[i-1] ];
        }


    printf("reading in %d lines ...\n",nrecs);

    /*try reading in one line at a time */
    for( j = 0, k = 0, l = 0; j < nrecs; j++) {
        for( i = 0; i < num; i++) {
            starts[i] = j;
            starts1[i] = l;
        }


        SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);

	/* get current position */
	x = *(double *)(btab + inoffsets[0]);
	y = *(double *)(btab + inoffsets[1]);
	z = *(double *)(btab + inoffsets[2]);
	vx = *(double *)(btab + inoffsets[3]);
	vy = *(double *)(btab + inoffsets[4]);
	vz = *(double *)(btab + inoffsets[5]);
	ident = *(int *)(btab + inoffsets[7]);


        SDFseekrdvecsarr(sdfp1, num, members, starts1, lines, addrs1, strides);

	/* get next position */
	x1 = *(double *)(btab1 + inoffsets[0]);
	y1 = *(double *)(btab1 + inoffsets[1]);
	z1 = *(double *)(btab1 + inoffsets[2]);
	vx1 = *(double *)(btab + inoffsets[3]);
	vy1 = *(double *)(btab + inoffsets[4]);
	vz1 = *(double *)(btab + inoffsets[5]);
	ident1 = *(int *)(btab1 + inoffsets[7]);

        dx = x1 - x;
        dy = y1 - y;
        dz = z1 - z;

        dvel = (vx1-vx)*dx + (vy1-vy)*dy + (vz1-vz)*dz;


        /* cut out 1. octant for more pleasant viewing experience */
        //if( (x < 0.0) || (y < 0.0) || (z < 0.0) ) {
        /* for all subsequent files, check if there are absorbed particles between 
           current and previous ident; pad with 0 D- and V- vector */

	if(first_file == 1) {
		idents[j] = *(int *)(btab + inoffsets[7]);
	}

        /* the already accreted particles */
	if(first_file == 0) {
            while(idents[k] < ident) {
               fprintf(fp,"D %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n",
                       0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f);
               fprintf(fp,"V %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n",
                       0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f);
	        ++k;
            }
        }

        /* put particles absorbed in sdfp1 at the origin; read in particles */
	while(ident < ident1) {
            /* the "offset vector". = x_n*(1+1/60) - x_(n+1)/60; for 60 fps */
            /* in this case, x+1 is 0,0,0 since particles are being accreted */
            offs_x =0.0;// x*1.01666666666666;
            offs_y =0.0;// y*1.01666666666666;
            offs_z =0.0;// z*1.01666666666666;

	    /* write initial position */
            fprintf(fp,"D %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n", 
                    x, -1.0*y, z, 0.f, 0.f, 1.f, 1.f);
            fprintf(fp,"V %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n",
                    offs_x, -1.0*offs_y, offs_z, 0.f, 0.f, 1.f, 1.f);

            ++j;
            ++k;
            for( i = 0; i < num; i++)
                starts[i] = j;
            SDFseekrdvecsarr(sdfp, num, members, starts, lines, addrs, strides);
            ident = *(int *)(btab + inoffsets[7]);
	    if( first_file == 1) 
                idents[j] = *(int *)(btab + inoffsets[7]);
            /* get next position */
            x = *(double *)(btab + inoffsets[0]);
            y = *(double *)(btab + inoffsets[1]);
            z = *(double *)(btab + inoffsets[2]);
	    vx = *(double *)(btab + inoffsets[3]);
	    vy = *(double *)(btab + inoffsets[4]);
	    vz = *(double *)(btab + inoffsets[5]);

            dx = x1 - x;
            dy = y1 - y;
            dz = z1 - z;
        
            dvel = (vx1-vx)*dx + (vy1-vy)*dy + (vz1-vz)*dz;

            /* in case any particles are accreted between ident and ident1 */
            if(first_file == 0) {
                while(idents[k] < ident) {
                    fprintf(fp,"D %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n", 
                            0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f);
                    fprintf(fp,"V %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n",
                            0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f);
    	            ++k;
                }
            }

	}

        offs_x = x1;//x*1.01666666666666 - x1*0.01666666666666;
        offs_y = y1;//y*1.01666666666666 - y1*0.01666666666666;
        offs_z = z1;//z*1.01666666666666 - z1*0.01666666666666;
        /* write current and next particle position */
        if(dvel > 0.0) {
            fprintf(fp,"D %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n", 
                    x, -1.0*y, z, 0.f, 0.f, 1.f, 1.f);
            fprintf(fp,"V %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n", 
                    offs_x, -1.0*offs_y, offs_z, 0.f, 0.f, 1.f, 1.f);
        } else {
            fprintf(fp,"D %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n", 
                    x, -1.0*y, z, 0.f, 0.f, 1.f, 1.f);
            fprintf(fp,"V %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f %+1.6f\n", 
                    offs_x, -1.0*offs_y, offs_z, 0.f, 0.f, 1.f, 1.f);
        }

        ++l;
        ++k;
        //}
    }

    printf("check: %d %d, %d %d, %d %d\n", nrecs, j, nrecs1, l, nrecs-nrecs1, k);

/*and we're done! clean up now -CE: if it ever works*/
/*    free(members);
    free(btab);
    free(addrs);
    free(types);
    free(inoffsets);
*/
    /*free(outbtab);*/
}


void writeheader(FILE *fp)
{
    // write output file header
    fprintf(fp, "set intensity EXPLICIT\n");
    fprintf(fp, "set parametric PARAMETRIC\n");
    fprintf(fp, "set filecontent DOTS\n");
    fprintf(fp, "set filetype NEW\n");
    fprintf(fp, "set depthcue 0\n");
    fprintf(fp, "set defaultdraw STELLAR\n");
    fprintf(fp, "set coordsys LEFT\n");
    fprintf(fp, "set author UNKNOWN\n");
    fprintf(fp, "set site UNKNOWN\n");
    fprintf(fp, "set library_id UNKNOWN\n");
    fprintf(fp, "set comment none\n");
}
