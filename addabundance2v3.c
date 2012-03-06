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
#include "cool.h"

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

//int NNW;

int calc_u;

//float **ionfracp;
//float **tablep;

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
    int **nwlist, Gridpts, Nel;
    char tmpchr[40], **names;
    char **vecs, **members, **outmembers;
    char **fnames, **pnames, **nnames;
    SDF_type_t *types, *outtypes;
    size_t stride = 0, outstride = 0;
    void *btab, *outbtab;
    void **addrs;
    float **abundarr; /*should this be void * ? */
    float *newabund, *radbin;
    float u_tot, n_e, eos_rho;
    double rho, temp, eos_n;
    double x, y, z, radius;

    init_CoolTable(&Gridpts, &Nel);
    printf("%d by %d cooling terms\n",Gridpts, Nel);

    frp = fopen("log.out","w");
    if(!frp) printf("error opening log file!\n");


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


        /* lastly, determine u from T. do in cgs */
        if(calc_u) {
            eos_rho = rho/(lengthCF*lengthCF*lengthCF)*massCF; /* cgs */
            eos_n = 0.0;
            for( i = 0; i < numA; i++) {
                eos_n += ((double)eos_rho)*N_AVOG /(double)(nparr[i] + nnarr[i]) *
                     (double)abundarr[jl][i];
            }
            n_e = find_ne(newabund, nparr, nnarr, temp, eos_rho, Gridpts, Nel); /* cgs */
            eos_n += n_e; /* code units */
            u_tot = (float)(1.5 * eos_n * K_BOLTZ * temp + 
                    A_RAD * temp * temp * temp * temp); 
            u_tot = (float)((double)u_tot/eos_rho); /* need specific internal energy */
            u_tot = u_tot * timeCF*timeCF/(lengthCF*lengthCF);
            memcpy(outbtab + inoffsets[ iu ], &u_tot, SDFtype_sizes[ types[ iu ]]);
        }


        /*dump the outbtab data into the file now-CE*/
        fwrite(outbtab, outstride, 1, fp);

    }
    printf("last: eos_n= %e n_e= %e\n", eos_n, n_e);
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

/* COPIED+PASTED FROM lcool1.c FROM THE SPH+NLN DIRECTORY */
/************************************************************************
 * PURPOSE:								*
 *   to calculate the cooling term for a given temperature and 		*
 *   composition, and return the cooling value to the calling routine.	*
 *   Returns total cooling term in erg*cm^3/s.				*
 *   The composition is taken directly from SPHbody->abund[i] and	* 
 *   converted from mass fraction to number density.			*
 *   									*
 * NOTE:								*
 *   expects table* and ionfrac* to have been initialized by 		*
 *   readintables.c							*
 *   the analytic cooling routine is still in here as analytic_cool	*
 *									*
 * DOES:								*
 *   - uses NR locate to find the index corresponding to the temperature*
 *   - extrapolate or return analytic cooling term if outside table	*
 *   - calculate cooling term for each element/ion and sum over those	*
 *   - uses NR polint to interpolate over the tables. 			*
 *   - it populates X_el from SPHbody->abund[i] and puts elements in 	*
 *     ascending order of Z, and also sums over isotopes. Thus the order*
 *     of isotopes/elements in SPHbody->abund does not matter. 		*
 ************************************************************************/

#include "physics_sph.h"
#include "cool.h"
#include "nrutil.h" /*ok to have this in*/

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif


/*
 * take all command line outputs like error or status messages out for 
 * now. eventually should have a flag in the .ctl files to turn 
 * debugging on (i.e. with output messages) or off
 */

/*arrays in C: array[row-index][col-index]*/


double find_ne(float abundarr[], int nparr[], int nnarr[] ,double temp, double rho, int Gridpts,int Nel)
{
/* could also calculate n_ion and set a flag that determines which one is returned */
    extern float **ionfracp;	/*global array with ion fractions*/
    double dy,df;	/*measure of accuracy returned from interp.*/
    double *dyp;	/*pointer to dy*/
    double *dfp;	/*pointer to df*/
    double fracn;	/*holds ion fraction returned by interp. */
    double *fracnp;	/*pointer to frac*/
    double fracneu;	/*fraction of neutral atoms per element*/
    double rowarr[2],fracns[2];	/*temporary arrays for interp.*/
    double logtemp; 	/*log(temp) for ionfracp interpolation*/
    float temps[Gridpts];
    float X_el[Nel];		/*element fraction, i.e. number density*/
    double nelectron = 0.;	/* total electron number density */
    long j;	/*holds index returned by locate routine*/
    long *jp;	/*pointer to j*/
    long j_prev;
    //int nparr[Nel], nnarr[Nel];	/*array for A,N numbers for each isotope*/
    int n,m,N,index; 	/*some indices for looping and arrays*/
	
    /*initialize things so I don't get stupid warnings all the time:*/
    j=-999;
    dy=-999;
    jp=&j;
    dyp=&dy;
    dfp=&df;
    fracnp=&fracn;

    /* we're doing temperature in logspace */
    logtemp=log10(temp);
    if( logtemp < ionfracp[0][0] ) logtemp = 4.0;
    if( logtemp > ionfracp[0][Gridpts-1] ) logtemp = 9.0;

    /*locate the indices of the table;
      same for both tables as they go over the same range/grid points 
    */
    for( n = 0; n < Gridpts; n++) {
        temps[n] = ionfracp[0][n];
    }

    locate(&temps[0], Gridpts, logtemp, jp); 
    j_prev=j;
	
    /*if we're outside the table, do analytic cooling if extrapolate=0
      or extrapolate if extrapolate=1:
      extrapolation is still a little funky - CE
      above table= j=-2, below table= j=-99
    */
    if (j==-2) j=0;
    if (j == -99) j=Gridpts-1; 

    for( n = 0; n < Nel; n++) X_el[n] = 0.; /*initialize all to zero*/

    for( n = 0; n < NISO; n++) {/*sum individual isotopes*/
        /* exclude bare neutrons;in tablep/ionfracp index=0 is H*/
	/* this also puts X_el in order of ascending Z, if an element
	 * does not exist in abund[i], it is just zero in X_el */
        if(nparr[n] > 0) {
           X_el[ nparr[n]-1 ] += abundarr[n] * rho *
               N_AVOG / ((float)(nnarr[n] + nparr[n]));
        }
    }

    for ( n = 0; n < Nel; n++)	
    {
	N = n+1;
	/*loop over ions for each element,dont skip any*/
        for ( m = 0; m < (N+1); m++) {
            index = (int)(N * (float)( (N+1)/2 ) ); /* to preempt int division probs */
	        /* re-assign table values so interpolation can be done in 
	           double precision (table is floats). */
	        rowarr[0] = ionfracp[0][ j ];
	        rowarr[1] = ionfracp[0][ j + 1 ];
		
	        fracns[0] = ionfracp[ index+m ][ j ];
	        fracns[1] = ionfracp[ index+m ][ j + 1 ];

	        /*interpolate*/
	        polint(rowarr,fracns,2,logtemp,fracnp,dfp);

	        /*reset value if extrapolated to unphysical value*/
	        if (fracn<0.0)
	        {
	            if ((fracns[0]-fracns[1]) <0) fracn=1.0e-40;
	    	    else fracn=1.0;
	        }

            nelectron += X_el[n] * (double)(m) * fracn;
	    }
    }
    return nelectron;
} /*end find_ne*/



double analytic_cool(double temp)
{
    /* From Chris's email; fit to Dalgarno and McCray (ARA&A
	 1972, 10, 375) and Sutherland and Dopita (ApJS, 88,
	 253) 
	 */
	double lcool;
	
	if (temp < 1.0e4)
		lcool = 1.0e-26 * exp(-1.0e5/temp) * sqrt(temp);//guessed; for O - CE
		//lcool = 1.0e-27 * exp(-1.0e2/temp) * sqrt(temp);
	else if (temp < 3.0e5)
		lcool=1.0e-21;
	else
		lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);
		//lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);
	
	return lcool;
}/*end analytic_cool*/



/* this bisection routine is from NR for C: */
/*"Given an array xx[1..n], and given an value x, returns a value
j such that x is between xx[j] and xx[j+1]. xx must be monotonic,
either increasing or decreasing. j=-2 or j=-99 is returned to indicate
that x is out of range.*/
/*#includes: none*/
/*call syntax: locate(&xx, N, x, j) */
/*~~~~~~~~~~~~~~~~WORKS!!! 04/23/2009~~~~~~~~~~~~~*/
void 
locate(float xx[], long Nel, float x, long *j)
{
   /*floats should be enough for finding correct indices*/
   long ju, jm, jl;
   int ascnd,sign;

   jl=-1;
   ju=Nel;
   /*check whether the table is in increasing (ascnd=1) or 
     decreasing (ascnd=0) order
   */
   ascnd=(xx[Nel-1] >= xx[0]);   /*what does this line do - check whether
                                   xx[N]>xx[1] or not -CE*/
   if (ascnd) sign=1;
   if (!ascnd) sign=-1;

   while (ju-jl > 1)
   {
      jm=(ju+jl) >> 1;
      if ((x >= xx[jm]) == ascnd)   
         jl=jm;
      else
         ju=jm;
   }
   /*if (true && true) && true (1 is true) then below (j<0) table*/
   if ( (((x - xx[0])*sign <0) && ((x-xx[Nel-1])*sign <0)) && 1)
	   *j=-2; 
   /*if (not false && not false) && true (1 is true) then above (j>Nel) table*/
   else if ( (!((x - xx[0])*sign <0) && !((x-xx[Nel-1])*sign <0)) && 1)
	   *j=-99; 
   else *j=jl;
} /*end locate*/



/*this is a polynomial interpolation routine from NR for C (S3.1):
"Given arrays xa[1..n] and and ya[1..n], and given a value x, this 
routine returns a value y, and an error estimate dy. If P(xP) is the 
polynomial of degree N-1 such that P(xa_i)=ya_i, i=1,...,n, then 
the returned value y=P(x).*/
/*#includes: math.h, "nrutil.h"*/
/*call syntax: polint(&xx[14],&yy[14],4,x,yp,dyp) for 4-point 
interpolation between tabulated points [14..17]*/
/*~~~~~~~~~~~~~~~~~ WORKS!!!! 04/23/2009 ~~~~~~~~~~~~~~~~~*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   //ZERO!!! offset is assumed in all indices
   int i,m,ns=0,size;
   double den,dif,dift,ho,hp,w;
   double *c, *d;

   //printf("interpolating between %E and %E .... get ", ya[0], ya[1]);
   //printf("xa= {%.1E, %.1E} \n", xa[0], xa[1]);
   //printf("ya= {%.1E, %.1E} \n", ya[0], ya[1]);

   //if we're exactly at one grid point, just return that value*/
   if ((x-xa[0]) == 0.0) 
	   *y=ya[0];
   else if ((x-xa[1]) == 0.0) 
	   *y=ya[1];
   else	//else interpolate*/
   {
   dif=fabs(x-xa[0]);
   c=dvector(0,n-1);
   d=dvector(0,n-1);
   //find the index ns of the closest table entry
   for (i=0;i<n;i++)
   {
      if ((dift=fabs(x-xa[i])) < dif)
      {
         ns=i;
         dif=dift;
      }
      //initialize the tableau of c's and d's
      c[i]=ya[i];
      d[i]=ya[i];
   }
   //initial approximation to y
   *y=ya[ns--];
   //for each column of the tableau ...
   for (m=0;m<n-1;m++)
   {
      //... loop over current c's and d's and update them
      for (i=0;i<n-m-1;i++)
      {
         ho=xa[i]-x;
         hp=xa[i+m+1]-x;
         w=c[i+1]-d[i];
	 den=ho-hp;
         //error message: two input xa's are identical to within roundoff
         //if (den==0.0) nerror("Error in routine polint");
         den=w/den;
         //update c's and d's
         d[i]=hp*den;
         c[i]=ho*den;
      }
      /*After each column in the tableau is completed, we decide wich 
        correction, c or d, we want to add to our accumulating value of
        y, i.e., which path to take through the tableau - forking up or 
        down. We do this in such a way as to take the most "straight 
        line" route through the tableau to its apex, updating ns 
        accordingly to keep track of where we are. This route keeps the 
        partial approximations centered (insofar as possible) on the 
        target x. the last dy added is thus the error indication. */
      *y += (*dy=(2*(ns+1)<(n-(m+1)) ? c[ns+1] : d[ns--]));
      //printf("%E\n",*y);
   }
   //printf("%E\n", *y);
   free_dvector(d,0,n-1);
   free_dvector(c,0,n-1);
   }
} //end polint



/*bilinear polynomial interpolation in 2 dimensions from NR for C: */
/*Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a 
submatrix of function values ya[1..m][1..n], tabulated at the grid points
defined by x1a and x2a; and given values x1 and x2 of the independent 
variables; this routine teturns as interpolated function value y, and an 
accuracy idication dy (based only on the interpolation in the x1 direction)*/
/*#includes: "nrutil.h" */
/* call syntax: polin2d(&x1a[jj],&x2a[kk],&yap[0][0],2,2,x1,x2,yp,dyp) 
 * where x1a,x1, refers to rows and jj is the row number from which to 
 * start the interpolation, x2a,x2 refers to columns and kk is the column
 * number from which to start the interpolation, 
 * and yap is an array of pointers (which point to the values
 * to interpolate) like so:
 * yap[0][0]=&table[jj][kk]
 * yap[1][0]=&table[jj][kk+1]
 * yap[0][1]=&table[jj+1][kk]
 * yap[1][1]=&table[jj+1][kk+1]
 * (really, whatever the order, polin2d treats the first index of yap as 
 * rows and interpolates over those first)
 */
/*~~~~~~~~~~~~~~~~ WORKS!!!! 04/23/2009 ~~~~~~~~~~~~~~*/
void
polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, 
        double x2, double *y, double *dy)
{
   //ZERO!!! offset is assumed in all indices
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

   int j;
   double *ymtmp;
 
   ymtmp=dvector(0,m-1);
   //loop over rows
   for (j=0;j<m;j++)
   {
      //interpolate over the 'rows' of the grid square containing the point
      //(x1,x2), i.e. over (x1a[j],x2). Put answer into temporary storage
      polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
   }
   //do the final interpolation, i.e. sort of interpolate over the 
   //interpolated 'rows'
   polint(x1a,ymtmp,m,x1,y,dy);
   free_dvector(ymtmp,0,m-1);
} //end polin2d


/* COPIED+PASTED FROM readintables.c FROM THE SPH+NLN DIRECTORY */
/******************************************************************************
 * the two arrays that hold the table values for the ion fractions (ionfracp) *
 * and cooling terms (tablep) are now 2D arrays, with the first index going   *
 * over the ions (order: H (Z=1) first, Zn (Z=30) last; for each element      *
 * neutral atom first, bare ion last), and the second index going over the    *
 * the data points for each temperature gridpoint. The correct row (ion) can  *
 * be found with: [ Nel*(Nel+1)/2 -1 + ionstate] ( -1 because arrays in C     *
 * start at 0), where ionstate is the ionization state (0 for neutral, Z+1    *
 * for bare ion), and Nel is the number of the element (=Z).                  *
 ******************************************************************************/
/*
extern float **tablep;
extern float **ionfracp;
*/

/* return Gridpts and Nel to calling function*/
void init_CoolTable(int *Gridpts, int *Nel)
{
    FILE *File1p;	/*pointer to file with cooling curves*/
    FILE *file2p;	/*pointer to file with ion fractions*/
    long lSize;		/*holds file size (number of characters in file)*/
    int i,j,k;		/*indices for looping through arrays*/
    int counter1 = 0, counter2 = 0;
    int index;		/*index to access correct array element*/
    int myint;		/*holds un-needed integers read in from files*/
    int tot_ion;	/* total number of ions in database (w/ bare ions)*/
    char mychar,myline[50];	/*holds new-lines and text read in from files*/

    /*open table with cooling curves*/
    File1p = fopen(
           "CHIANTI-COOLING.dat", "r");
           /*"/home/cellinge/SNSPH.dir/tree16/sph+nln/CHIANTI-COOLING.dat", "r");*/
    if (File1p == NULL) 
	singlPrintf("error opening cooling curves: \n");

    /*open table with ion fractions*/
    file2p = fopen(
           "mazzotta_etal_9.ioneq","r");
           /*"/home/cellinge/SNSPH.dir/tree16/sph+nln/mazzotta_etal_9.ioneq","r");*/
    if (file2p == NULL) 
	singlPrintf("error opening ion fractions: \n");

    fscanf(file2p, "%i %i", Gridpts,Nel);
     
    /* (Nel*(Nel+1)/2 + Nel is total number of ions, incl. bare ions */
    tot_ion = (*Nel) * ( (*Nel)+3 ) / 2;
    /*singlPrintf("total ions: %d\n", tot_ion);*/


    /* for ionfracp, need (Nel*(Nel+1)/2 + Nel+1) by Gridpts array */
    ionfracp = (float**)malloc( (tot_ion + 1) * sizeof(float *) );

    for ( i = 0; i < (tot_ion + 1); i++) {
	ionfracp[i] = (float *)malloc( (*Gridpts) * sizeof(float) );
        for ( j = 0; j < (*Gridpts); j++)
            ionfracp[i][j] = 0.0;
    }


    /* for tablep, need (Nel*(Nel+1)/2) by Gridpts array */ 
    tablep = (float **)malloc( (tot_ion) * sizeof(float *) );

    for ( i = 0; i < (tot_ion); i++) {
	tablep[i] = (float *)malloc( (*Gridpts) * sizeof(float *) );
        for ( j = 0; j < (*Gridpts); j++)
            tablep[i][j] = 0.0;
    }

    fgets(myline,50,File1p);//read in first line of text in cooling curves
    fgets(myline,50,File1p);//read in second line of text in cooling curves

    /* since the ions are from H to Zn in ascending order, don't need
     * Z of element (?) */
    for( i = 0; i < (*Nel); i++) 
	    fscanf(File1p, "%*i");


    mychar=fgetc(File1p);//read in extra new-line in cooling curves
    fgets(myline,50,File1p);//read in line "temperatures...." in cooling curves


    /* we're getting log(T) from ionfractions, so skip T from cooling table*/
    for ( i = 0; i < (*Gridpts); i++) 
	    fscanf(File1p, "%*g");

    /*read in log(temperatures) for ion fractions*/
    for ( i = 0; i < (*Gridpts); i++) 
	    fscanf(file2p,"%4g",&ionfracp[0][i]);
	

    fgets(myline,50,File1p);//get trailing new-line in cooling curves
    fgets(myline,50,File1p);//get trailing new-line in cooling curves


    /*loop over elements*/
    for ( j = 0; j < (*Nel); j++) { 
	/*get line with element name*/
        fgets(myline,50,File1p);

        /*loop over ions in current element*/
	for ( k = 0; k <= (j+1); k++) {
	    fscanf(file2p,"%3i",&myint);
	    fscanf(file2p,"%3i",&myint);

            for ( i = 0; i < (*Gridpts); i++) {

    	        fscanf(File1p,"%13E",&tablep[ counter1 ][i]);

		/*add one to counter2 since the 1st contains log(T) */
	        fscanf(file2p,"%10E",&ionfracp[ counter2+1 ][i]);

	/* in CHIANTI-file, (neutrals?) bare ions seem to be missing, but are
	 * present in Mazzotta et al. ion fractions file. In the above set-up
	 * the bare ions from the ion-fractions file are read in, and the ions
	 * in tablep and ionfracp are lined up (i.e. each row of ionfrac
	 * and tablep contain the same ion of the same element), and the 
	 * bare ions are read in as zeros in tablep).
	 */
	    }

            //printf("counter: %4d, j=%2d, k=%2d, tab= %.6E  ion=%.6E\n",counter1,j,k,tablep[counter1][12],ionfracp[counter1][12]);
	    counter1++;
	    counter2++;
	}
	mychar=fgetc(File1p);//get first trailing new-line
	mychar=fgetc(File1p);//get second trailing new-line
    }

    //close file
    fclose(File1p);
    fclose(file2p);
    /*DO NOT free(tablep); UNTIL THE END OF THE WHOLE PROGRAM!!!!*/
}

