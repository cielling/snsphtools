/*
   PURPOSE:
        to make a table correlating particle id with zone from abundance table.

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

static const int NISO = 22;

static void initargs(int argc, char *argv[], SDF **sdfp)
{
    *sdfp = SDFopen(NULL, argv[1]);
    if (*sdfp == NULL) {
	fprintf(stderr, "%s: %s: %s\n", argv[0], argv[1], SDFerrstring);
	exit(2);
    }
}

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;
    FILE *f2p = NULL;
    int nmembers, nrecs, nvecs, idindex, ident, i, j, k, jl, jm, ju;
    int numA, Nbins, len, counter, flag = 0;
    int *nparr, *nnarr;
    char mychar[5], **names;
    size_t stride = 0, outstride = 0;
    float **abundarr; /*should this be void * ? */
    double x, y, z, radius;
    float *radbin;
    char **vecs, **members;
    int *strides, *nobjs, *starts, *inoffsets;
    void **addrs, *btab;
    SDF_type_t *types;

    initargs(argc, argv, &sdfp);

    SDFgetint(sdfp, "npart", &nrecs);
    printf("%d particles\n",nrecs);

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) ++nmembers;
    }

/*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc( nmembers * sizeof (char *));
    addrs = (void **)malloc( nmembers * sizeof( void * ));
    strides = (int *)malloc( nmembers * sizeof(int));
    nobjs = (int *)malloc( nmembers * sizeof(int));
    starts = (int *)malloc( nmembers * sizeof(int));
    types = (SDF_type_t *)malloc( nmembers * sizeof( SDF_type_t ));
    inoffsets = (int *)malloc( nmembers * sizeof(int));

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
        if (strncmp(vecs[i], "ident", strlen(vecs[i])) == 0) idindex = nmembers-1;
    }

    btab = (void *)malloc(stride);

    /*calculate the byte offset in memory to the address of the next member-CE*/
    for (i = 0; i < nmembers; ++i) {
	addrs[i] = (char *)btab + inoffsets[i];
	strides[i] = stride;
    }


    /* now need to read in the abundances, and other necessary stuff */

    fp = fopen("abund.txt", "r");
    if (!fp) printf("error opening file abund.txt\n");

    /*read in first number, that's the number of isotopes in the file*/
    fscanf(fp, "%d", &numA);

    /*read in second number, that's the number of lines in the file*/
    fscanf(fp, "%d", &Nbins);

    nparr = (int *)malloc(numA * sizeof(int));
    nnarr = (int *)malloc(numA * sizeof(int));

    for ( i = 0; i < numA; i++) fscanf(fp, "%d", &nparr[i]);
    for ( i = 0; i < numA; i++) fscanf(fp, "%d", &nnarr[i]);

    names = (char **)malloc((numA+1) * sizeof(char *) );
    for ( i = 0; i < numA + 1; i++) {
        names[i] = (char *)malloc(4*sizeof(char));
        fscanf(fp, "%s", &names[i]);
    }

    /* malloc 2D array for abundances[radial bin][isotope] -CIE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    /* read in radial bins into radbin */
    radbin = (float *)malloc(Nbins * sizeof(float));
    if(!radbin) printf("error allocating radbin\n");

    /*read in abundances, one set for each radial bin*/
    for (i=0; i< Nbins; i++){
        fscanf(fp,"%12G", &radbin[i]);
        for (j=0; j < numA; j++) {
            fscanf(fp, "%12G", &abundarr[i][j]);
        }
    }

    fclose(fp);
    fp = NULL;


    fp = fopen("zonetopartid.dat", "w");
    if(!fp) printf("error opening file 'zonetopartid.dat'\n");

    for (j = 0; j < nrecs; ++j) {
        //printf("%8d ",j);
        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs, strides);

        /* calculate radius. this should work -CE */
        x = *((double *)(btab + inoffsets[0]));
        y = *((double *)(btab + inoffsets[1]));
        z = *((double *)(btab + inoffsets[2]));
        radius = sqrt( x*x + y*y + z*z );

        ident = *((int *)(btab + inoffsets[ idindex ]));

        /*quick bisection to locate the radial bin I'm in -CE */
        /* put everything below radbin[0] into first bin (bin zero),
           and put everything above radbin[Nbins-1] into last bin */
        if (radius >= radbin[Nbins-1])
            jl = Nbins-1;
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

        //printf("rad= %E id= %d zone= %d\n", radius, ident, jl);

/* which index holds my bin number? -CE: jl */

        fprintf(fp,"%10d ", ident);
        fprintf(fp,"%10d ",jl);
        fprintf(fp,"\n");

	for (i = 0; i < nmembers; ++i) {
            starts[i] = j+1;
	}
    }

    fclose(fp);

    return 0;
}
