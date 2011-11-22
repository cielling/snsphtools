/*

   PURPOSE:
	to add an abundance table to an initial conditions sdf file.
	it is assumed that fewer isotopes are added than are in the abundance
	files, so the added abundances are renormalized to preserve the
	original neutron-proton ratio.
    The neutron excess is stuffed into 56Fe, and the rest of the abundances
    is adjusted accordingly so they still sum to one.
        This version should also work for large files, since it reads in the
        sdf file line-by-line.
	New: this version writes the new abundance format, with the Z,N of the
	tracked isotopes in the header of the SDF file. It also expects to find
    a file "network.isotopes" with the Z,N numbers of the isotopes in the
    network in tabular form, and the first line containing the number of
    isotopes.

   COMPILE:
	with Makefile2. un-comment appropriate lines.
	This routine needs some libraries from the tree code and SDF routines,
	so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
	once for serial use (i.e. without PAROS flag).

   RUN:
	addabundance2new <in-file.sdf> <out-file.sdf>

   METHOD:
	@ read in <in-file.sdf> line-by-line (particle-by-particle)
	@ read in abundance information file
	@ write new header with added columns (struct members) for added
	  abundances to <out-file.sdf>
	@ calculated radius of each particle read in
	@ find correct radial bin
	@ renormalize abundances to retain approximate neutron excess
	@ write particle with selected abundance to <out-file.sdf>

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <SDF.h>
#include "consts.h"


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

int NNW;
int NNW;

int calc_u;

/*void renorm(int want[][NNW], float newabund[], float abundarr[], int nparr[], int nnarr[], int Niso, int Nnwi, FILE *frp);*/
void renorm(int **want, float newabund[], float abundarr[], int nparr[], int nnarr[], int Niso, int Nnwi, FILE *frp);

double get_ye(float abundarr[], int nparr[], int nnarr[], int Niso);

int make_spec_names(char *** chararr, char spec, int num);

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

/* open/initialize all files */
static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp)
{
    char input;

    if (argc != 4) {
	fprintf(stderr, "Usage: %s SDFfile outfile calc_u_flag\n", argv[0]);
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

    calc_u = atoi(argv[3]);
    if(calc_u) printf("calculating u from T\n");
    else printf("leaving u unchanged\n");
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

    /*figure out number of lines in the header, basically-CE*/
    nvecs = SDFnvecs(sdfp);
    /*get the names of the variables/parameters in the header-CE*/
    vecs = SDFvecnames(sdfp);

    /*go through all of them individually-CE*/
    for (i = 0; i < nvecs; ++i) {

        /* figure out where the data columns start - don't think this is necessary yet */
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
    FILE *afp = NULL, *frp = NULL;
    int i, j, k, nvecs, nmembers, Amembers, ju, jm, jl;
    int numA, Nbins, counter, flag = 0;
    int nrecs = 0, irho, iu, itemp;
    int *strides, *nobjs, *starts, *inoffsets, *outoffsets, *nparr, *nnarr;
    int **nwlist;
    char tmpchr[40], **names;
    char **vecs, **members, **outmembers;
    char **fnames, **pnames, **nnames;
    SDF_type_t *types, *outtypes;
    size_t stride = 0, outstride = 0;
    void *btab, *outbtab;
    void **addrs;
    float **abundarr; /*should this be void * ? */
    float *newabund, *radbin;
    float u_tot;
    double kb, arad, rho, temp, eos_n;
    double x, y, z, radius;

    frp = fopen("log.out","w");
    if(!frp) printf("error opening log file!\n");

    kb=K_BOLTZ *(timeCF*timeCF)/(massCF *lengthCF*lengthCF);
    arad=A_COEFF*(lengthCF*timeCF*timeCF/massCF);

    SDFgetint(sdfp, "npart", &nrecs);
    printf("%d particles\n", nrecs);

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    /* this assumes that "x" is the first data column. can't use nrecs,
       since that read in the whole file, which we're trying to avoid */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
		if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

    /* read in the list of isotopes in the network */
    afp = fopen("network.isotopes", "r");
    if (!afp) printf("error opening files network.isotopes\n");

    fscanf(afp, "%d", &NNW);
    fscanf(afp, "%s\t%s", tmpchr,tmpchr);

    /* this way nwlist is indexed the same way as want */
    nwlist = (int **) malloc ( 2 * sizeof(int *) );
    nwlist[0] = (int *) malloc ( NNW * sizeof(int) );
    nwlist[1] = (int *) malloc ( NNW * sizeof(int) );

    for( i = 0; i < NNW; i++) {
        fscanf(afp, "%d", &nwlist[0][i]);
        fscanf(afp, "%d", &nwlist[1][i]);
    }

    fclose(afp);
    afp = NULL;

    //NNW = NNW;
    Amembers = 1*NNW;

    /*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    outmembers = (char **)malloc(1 * NNW * sizeof(char *));
    for(j=0;j<1*NNW;j++) outmembers[j] = (char *)malloc(10 * sizeof(char));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    strides = (int *)malloc(nmembers * sizeof(int));
    nobjs = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    outtypes = (SDF_type_t *)malloc(1 * NNW * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));
    outoffsets = (int *)malloc(Amembers * sizeof(int));

    /*one by one, go through all the fields in the column, i.e. members of the struct?-CE*/
    flag = 0;
    for (i = 0, nmembers = 0, stride = 0; i < nvecs; ++i) {
        if( strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) irho=i;
        if( strncmp(vecs[i], "u", strlen(vecs[i])) == 0) iu=i;
        if( strncmp(vecs[i], "temp", strlen(vecs[i])) == 0) itemp=i;
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	    if (flag) {
	        members[nmembers] = vecs[i];
	        nobjs[nmembers] = 1;/*nobjs[0] is the number of particles. are all elements
                                      in this nobjs array the same then??-CE:yes, but can be
                                      different*/
			/* this is also the number of lines that are read in at a time */
	        starts[nmembers] = 0;  /* Not correct in parallel; use
				      NobjInitial */
	        types[nmembers] = SDFtype(members[nmembers], sdfp);
	        inoffsets[nmembers] = stride;
	        stride += SDFtype_sizes[types[nmembers]];

	        ++nmembers;
	    }
    }

    irho -= (nvecs - nmembers);
    iu -= (nvecs - nmembers);
    itemp -= (nvecs - nmembers);

    btab = (void *)malloc(stride);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
    for (i = 0; i < nmembers; ++i) {
	addrs[i] = (char *)btab + inoffsets[i];
	strides[i] = stride;
    }

/* up to here just the SDF file is read in, so everything stays the same -CIE */

/*
    fnames = make_spec_names('f', NNW);
    pnames = make_spec_names('p', NNW);
    nnames = make_spec_names('n', NNW);
*/
    /* what to call mass fraction, Z and N of isotope */
    make_spec_names(&fnames, 'f', NNW);
    make_spec_names(&pnames, 'p', NNW);
    make_spec_names(&nnames, 'n', NNW);

    outstride = stride;

    /*calculate at what byte-intervals the data should be written-CIE*/
    /* adding enough columns for the abundances */
    for (i = 0; i < Amembers; ++i) {
        sprintf(outmembers[i],"%s",fnames[i]);
        outtypes[i] = SDF_FLOAT;
	    outoffsets[i] = outstride;
	    outstride += SDFtype_sizes[ outtypes[ i ] ];
    }

    /*malloc enough space in memory for the array that holds the whole output data-CIE*/
    printf("about to malloc %d bytes for outbtab\n", (int)outstride);
    outbtab = (void *)malloc( outstride );

    /* now need to read in the abundances, and other necessary stuff */

    /* get Z and N of isotopes - WORKS!!*/
/*
    afp = fopen("c16r4_abun.dat", "r");
*/
    afp = fopen("abund.txt", "r");
    if (!afp) printf("error opening file abund.txt\n");

    /*read in first number, that's the number of isotopes in the file*/
    fscanf(afp, "%d", &numA);

    /*read in second number, that's the number of lines in the file*/
    fscanf(afp, "%d", &Nbins);

    nparr = (int *)malloc(numA * sizeof(int));
    nnarr = (int *)malloc(numA * sizeof(int));

    for ( i = 0; i < numA; i++) fscanf(afp, "%d", &nparr[i]);
    for ( i = 0; i < numA; i++) fscanf(afp, "%d", &nnarr[i]);

    names = (char **)malloc((numA+1) * sizeof(char *) );
    for ( i = 0; i < numA + 1; i++) {
        names[i] = (char *)malloc(4*sizeof(char));
        fscanf(afp, "%s", &names[i]);
    }

    /* malloc 2D array for abundances[radial bin][isotope] -CIE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    newabund = (float *)malloc(NNW * sizeof( float ));
    if (!newabund) printf("error allocating newabund\n");

    /* read in radial bins into radbin */
    radbin = (float *)malloc(Nbins * sizeof(float));
    if(!radbin) printf("error allocating radbin\n");

    /*read in abundances, one set for each radial bin*/
    for (i=0; i< Nbins; i++){
        fscanf(afp,"%12G", &radbin[i]);
        for (j=0; j < numA; j++) {
            fscanf(afp, "%12G", &abundarr[i][j]);
        }
    }

    fclose(afp);
    afp = NULL;


    /* print Z and N of isotopes of choice */
    for( i = 0; i < NNW; i++ ) {
        fprintf(fp, "int %s = %d;\n", pnames[i], nwlist[0][i]);
        fprintf(fp, "int %s = %d;\n", nnames[i], nwlist[1][i]);
    }

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
    fprintf(fp, "}[%d];\n", nrecs);
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
        rho = *(float *)(btab + inoffsets[irho]);
        temp = *(float *)(btab + inoffsets[itemp]);

        /*quick bisection to locate the radial bin I'm in -CE */
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
        /* which index holds my bin number? -CE: jl */

        /* copy data from sdf file to output array */
	for (i = 0; i < nmembers; ++i) {
            memcpy(outbtab + inoffsets[i], btab + inoffsets[i],
	           SDFtype_sizes[ types[i] ]);
            starts[i] = j+1;
	}

        /* renormalize abundances first */
        renorm(nwlist, newabund, abundarr[jl], nparr, nnarr, numA, NNW, frp);

        /* now populate outbtab with abundances 'n stuff */
        for (i = 0, counter = 0; i < numA; i++) {
            for (k = 0; k < NNW; k++) {
                if ( (nparr[i] == nwlist[0][k]) &&
                     (nnarr[i] == nwlist[1][k]) ) {

                     /* fill in abundances */
                    memcpy(outbtab + outoffsets[ k ],
                       &newabund[k] , SDFtype_sizes[ outtypes[ k ] ]);
                       /*&abundarr[jl][i] , SDFtype_sizes[ outtypes[ k ] ]);*/

                    counter++;
                }
            }
        }


        /* lastly, determine u from T */
        if(calc_u) {
            eos_n = 0.0;
            for( i = 0; i < numA; i++) {
                eos_n += ((double)rho)*N_AVOG * massCF /(double)(nparr[i] + nnarr[i]) *
                     (double)abundarr[jl][i];
                     // * (double)(nparr[j] + 1.0);/* accounts for electrons!*/
            }
            u_tot = (float)(1.5 * eos_n * kb * temp + 
                    arad * temp * temp * temp * temp); 
            u_tot = (float)((double)u_tot/rho); /* need specific internal energy */
            memcpy(outbtab + inoffsets[ iu ], &u_tot, SDFtype_sizes[ types[ iu ]]);
        }


        /*dump the outbtab data into the file now-CE*/
        fwrite(outbtab, outstride, 1, fp);

    }
    printf("wrote %d abundances to outbtab\n", counter);

    /*and we're done! clean up now -CE*/
    free(pnames);
    free(nnames);
    free(fnames);
    printf("free'd names\n");
    free(members);
    free(addrs);
    free(strides);
    printf("free'd members addrs strides\n");
    free(nobjs);
    free(starts);
    free(types);
    printf("free'd nobjs starts types\n");
    free(inoffsets);
    free(outoffsets);
    printf("free'd offsets\n");

    free(btab);
    free(outbtab);
    printf("free'd btabs\n");
}


double get_ye(float abundarr[], int nparr[], int nnarr[], int Niso) {
    int i;
    double Ye;
    Ye = 0.0;
    for( i=0; i<Niso; i++) {
        Ye += (double)abundarr[i] / (double)( nparr[i]+nnarr[i] ) * (double)nparr[i];
    }
    return Ye;
}


/*void renorm(int want[][NNW], float newabund[], float abundarr[], int nparr[], int nnarr[], int Niso, int Nnw, FILE *frp) */
void renorm(int **want, float newabund[], float abundarr[], int nparr[], int nnarr[], int Niso, int Nnw, FILE *frp)
{
    int i, j, is_in_arr, adjust;
    int jl, jm, ju, fe56, count;
    double Ye_old, Ye_new, sum, factor, factor1, eps=1.e-4;

    for( i=0; i<Nnw; i++) newabund[i] = 0.0;

    for( j=0; j<Niso; j++) {

        is_in_arr = 0;

        for( i=0; i<Nnw; i++) {
            if((want[0][i] == nparr[j]) && (want[1][i] == nnarr[j])) {
                newabund[i] += abundarr[j];
                is_in_arr = 1;
            }
            if((want[0][i] == 26) && (want[1][i] == 30)) fe56 = i;
        }

		/* figure out where we are */
        if(is_in_arr == 0) {
            if (nparr[j] >= want[0][Nnw-1])
                jl = Nnw-1;
            else if (nparr[j] <= want[0][0])
                jl = 0;
            else {
                ju = Nnw - 1;
                jl = 0;
                while (ju-jl > 1) {
                   jm = (ju+jl) >> 1;
                   if ((nparr[j] >= want[0][jm]) == (want[0][Nnw-1] >= want[0][0]))
                      jl=jm;
                   else
                      ju=jm;
                }
            }

            /* determine appropriate abundance bin to stuff isotope into */
            if(jl == 0) {
                newabund[jl] += abundarr[j];
            } else if(jl == Nnw-1) {
                newabund[jl] += abundarr[j];
            } else {
                if(want[0][jl] >= nparr[j]) {
                    if((want[0][jl]-nparr[j]) <= (nparr[j]-want[0][jl-1])) {
                        newabund[jl] += abundarr[j];
                        //printf("%2d: p=%2d\n",jl,want[0][jl]);
                    } else {
                        newabund[jl-1] += abundarr[j];
                        //printf("%2d: p=%2d\n",jl-1,want[0][jl-1]);
                    }
                } else {
                    if((want[0][jl+1]-nparr[j]) <= (nparr[j]-want[0][jl])) {
                        newabund[jl+1] += abundarr[j];
                        //printf("%2d: p=%2d\n",jl+1,want[0][jl+1]);
                    } else {
                        newabund[jl] += abundarr[j];
                        //printf("%2d: p=%2d\n",jl,want[0][jl]);
                    }
                }
            }
        }
    }

    sum = 0.0;
    for( i=0; i<Nnw; i++) sum += newabund[i];
    if(fabs(sum - 1.0) > 1.e-3)
       fprintf(frp, "warning: mass fractions don't sum to 1: %E\n",sum);

    /* check and adjust the Ye's */
    count = 0;
    do {
        Ye_old = get_ye(abundarr, nparr, nnarr, Niso);
        Ye_new = get_ye(newabund, want[0], want[1], Nnw);

        if(fabs(Ye_old/Ye_new - 1.0) > eps) {

            /* new Ye too small: decrease Fe56. new Ye too large: increase Fe56 */
            factor = Ye_new/Ye_old;
            factor1 = ( (1.0 - factor*newabund[fe56]) / (1.0-newabund[fe56]));
            newabund[fe56] *= factor;

            for( i=0; i<Nnw; i++) {
                if(i != fe56) {
                   newabund[i] *= factor1;
                   if(newabund[i] < 0.000) fprintf(frp,
                      "warning: negative abundance for %2d %2d after %d iterations!\n",
                      want[0][i],want[1][i],count);
                }
            }

            adjust = 1;
            count++;

        } else {
            adjust = 0;
        }

        if(count > (int)(2.0/eps)) adjust = 0;
            adjust = 0;

    } while (adjust == 1);

    fprintf(frp, "adjusted Ye from %E to %E\n",Ye_old, Ye_new);

    sum = 0.0;

    for( i=0; i<Nnw; i++) sum += newabund[i];

    if(fabs(sum - 1.0) > 1.e-3)
       fprintf(frp, "renorm warning: mass fractions don't sum to 1: %E\n",sum);

}


int make_spec_names(char ***chararr, char spec, int num)
{
    int i;
    char tmpchr[20];

    *chararr = (char **)malloc(num * sizeof(char *) );

    for( i = 0; i < num; i++ ){

        sprintf( tmpchr, "%c%-d\0", spec, (i+1));
        (*chararr)[i] = (char *)malloc( ( strlen(tmpchr) + 1) * sizeof(char) );
        sprintf( (*chararr)[i], "%s",tmpchr );
        //printf("isotope specifier: %s  %d\n", specarr[i],strlen(specarr[i]));

    }
    return i;
}
