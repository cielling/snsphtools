/*
   PURPOSE:
	to add a velocity asymmetry to an sdf file
        This version should also work for large files, since it reads in the
        sdf file line-by-line.

   COMPILE:
    with Makefile2. un-comment appropriate lines.
    This routine needs some libraries from the tree code and SDF routines,
    so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
    once for serial use (i.e. without PAROS flag).

   RUN:
    addvel <in-file.sdf> <out-file.sdf>

   METHOD:
	@ read in header of <in-file.sdf>
	@ prompt user for the shock radius (or below where-ever the velocity
	  asymmetry is to be added
	@ read in data line-by-line
	@ calculate radius for each particle
	@ determine if particle is in/inside of the shock (by simple comparison
      to the input radius)
	@ if yes, modify velocity according to prescription
	@ write data to <out-file.sdf>

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
static void writeinit(FILE *fp);
static void writescalars(SDF *sdfp, FILE *fp);
static void writestructs(SDF *sdfp, FILE *fp);

static const int NISO = 22;
static int asym_u = 0;

void renorm(int want[][NISO], float newabund[], float abundarr[], int nparr[], int nnarr[], int Niso, int Nnwi, FILE *frp);

double get_ye(float abundarr[], int nparr[], int nnarr[], int Niso);

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

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp)
{
    char input;

    if (argc != 4) {
	fprintf(stderr, "Usage: %s SDFfile outfile asym_u\n", argv[0]);
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

    asym_u = atoi(argv[3]);
    if(asym_u == 1) printf("calculating asymmetry in u also\n");
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

    nvecs = SDFnvecs(sdfp);/*figure out number of lines in the header, basically-CIE*/
    vecs = SDFvecnames(sdfp);/*get the names of the variables/parameters in the header-CIE*/

    /*go through all of them individually-CIE*/
    for (i = 0; i < nvecs; ++i) {
        /*figure out if scalar(=1) or array(!=1)?-CIE*/
/*
	nrecs = SDFnrecs(vecs[i], sdfp);
	if (nrecs != 1) continue;
*/
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) continue;

	type = SDFtype(vecs[i], sdfp);
	/* Doesn't handle string */
	/* Also, SDFgetfloat was happy to read a double scalar and
	   convert it for me; that probably works backwards too.  I
	   don't think there's an equivalent for the SDFrdvecs family
	   though, so "read; convert type; write" is the general
	   path. */
        /*read in header file, line by line, with the appropriate function-CIE*/
	switch (type) {
	case SDF_INT:
        /*SDFget*:sdfp=sdf file; vecs[i]= variable name; datum.*=holds value of that variable-CIE*/
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

        /*write header file, line by line, as the appropriate data type-CIE*/
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
    FILE *frp = NULL;
    int i, j, nvecs, nmembers;
    int len, flag = 0;
    int nrecs;
    int ixvel, iu, ih;
    int *strides, *nobjs, *starts, *inoffsets;
    char **vecs, **members;
    SDF_type_t *types;
    size_t stride = 0;
    void *btab;
    void **addrs;
    double alpha, beta, vfactor, cos_angle;
    float set_radius;
    double x, y, z, radius;
    float vx, vy, vz, vx2, vy2, vz2, v_r, u;

    frp = fopen("/work/01834/cielling/log.out", "w");
    if(!frp) printf("error opening log file!\n");

    /* j2 geometry; from Hungerford, Fryer, & Warren 2003 */
/*
    alpha = sqrt(3./7.);
    beta = sqrt(3./7.);
*/

    /* j4 geometry */
/*
    alpha = sqrt(2./7.)-0.2;
    beta = sqrt(9./7.);
*/

    /* e2 geometry */
/*
    alpha = 4./3.;
    beta = -2./3.;
*/

    /* e4 geometry */
    /* subtrahend is an attempt to conserve energy better */
    alpha = 8./5. -0.05;
    beta = -6./5.;

    printf("set_radius: enter 0 if set by vr:");
    scanf("%f",&set_radius);
    printf("%f\n",set_radius);

    vfactor = 1.00;

    SDFgetint(sdfp, "npart", &nrecs);
    fprintf(stdout,"%d particles\n", nrecs);

    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    /* Count structure members */
    for (i = 0, nmembers = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) ++nmembers;
    }
    fprintf(stdout,"nmembers = %d\n",nmembers);


    /*malloc memory space for the respective features of the struct-CIE*/
    members = (char **)malloc(nmembers * sizeof(char *));
    addrs = (void **)malloc(nmembers * sizeof(void *));
    strides = (int *)malloc(nmembers * sizeof(int));
    nobjs = (int *)malloc(nmembers * sizeof(int));
    starts = (int *)malloc(nmembers * sizeof(int));
    types = (SDF_type_t *)malloc(nmembers * sizeof(SDF_type_t));
    inoffsets = (int *)malloc(nmembers * sizeof(int));

    /*one by one, go through all the fields in the column, i.e. members of the struct?-CIE*/
    flag = 0;
    for (i = 0, nmembers = 0, stride = 0; i < nvecs; ++i) {
        if (strncmp(vecs[i], "x", strlen(vecs[i])) == 0) flag = 1;
	if (flag) {
	    members[nmembers] = vecs[i];
	    nobjs[nmembers] = 1;/*nobjs[0] is the number of particles. are all elements
                                      in this nobjs array the same then??-CIE:yes, but can be
                                      different*/
	    starts[nmembers] = 0;  /* Not correct in parallel; use
				      NobjInitial */
	    types[nmembers] = SDFtype(members[nmembers], sdfp);
	    inoffsets[nmembers] = stride;
            if (strncmp(vecs[i], "vx", strlen(vecs[i])) == 0) ixvel = nmembers;
            if (strncmp(vecs[i], "u", strlen(vecs[i])) == 0) iu = nmembers;
            if (strncmp(vecs[i], "h", strlen(vecs[i])) == 0) ih = nmembers;
	    stride += SDFtype_sizes[types[nmembers]];

	    ++nmembers;
	}
    }

    fprintf(stdout,"vx at %d\n",ixvel);

    btab = (void *)malloc(stride);

    /*calculate the byte offset in memory to the address of the next member-CIE*/
    for (i = 0; i < nmembers; ++i) {
	addrs[i] = (char *)btab + inoffsets[i];
	strides[i] = stride;
    }


    /*print the struct declaration part from the header-CIE*/
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
    fprintf(fp, "}[%d];\n", nrecs);
    fprintf(fp, "#\n");
    fprintf(fp, "# SDF-EOH\n");

    for (j = 0; j < nrecs; ++j) {

        for( i = 0; i < nmembers; i++) starts[i] = j;
        //printf("%8d ",j);
        /* read one line of data into btab (via addrs) */
        SDFseekrdvecsarr(sdfp, nmembers, members, starts, nobjs, addrs,
		     strides);

        /* calculate radius. this should work -CIE */
        x = *(double *)(btab + inoffsets[0]);
        y = *(double *)(btab + inoffsets[1]);
        z = *(double *)(btab + inoffsets[2]);
        radius = sqrt( x*x + y*y + z*z );
        if( beta < 0.) 
            cos_angle = fabs(x)/radius; /* equatorial asymmetry */
        else 
            cos_angle = fabs(z)/radius; /* polar (jet) asymmetry */

        /* get velocities */
        vx = *(float *)(btab + inoffsets[ixvel]);
        vy = *(float *)(btab + inoffsets[ixvel+1]);
        vz = *(float *)(btab + inoffsets[ixvel+2]);
        v_r = ((vx*x)+(vy*y)+(vz*z))/radius;

        /* get internal energy */
        if(asym_u == 1) {
            u = *(float *)(btab + inoffsets[iu]);
            h = *(float *)(btab + inoffsets[ih]);
            h3 = h*h*h;
        }

        theta = (float)(fabs(z*z)/(radius*radius)); /* negative for eq. asymmetry */

        /* calculate the velocity asymmetry */
        /* only for expanding velocities */
        if( (v_r > 0.0 && set_radius == 0)
            || (radius < set_radius)) {
            /* jet geometry */
            vx2 = (alpha + beta * cos_angle) * vx;
            vy2 = (alpha + beta * cos_angle) * vy;
            vz2 = (alpha + beta * cos_angle) * vz;
            if(asym_u == 1) {
                u = (alpha + beta * cos_angle) * sqrt(u);
                u = u*u; /* to conserve energy, since k.e. ~ v^2 and using same formula */
                memcpy(btab + inoffsets[iu], &u, SDFtype_sizes[ types[iu] ]);
            }

            memcpy(btab + inoffsets[ixvel], &vx2, SDFtype_sizes[ types[ixvel] ]);
            memcpy(btab + inoffsets[ixvel+1], &vy2, SDFtype_sizes[ types[ixvel+1] ]);
            memcpy(btab + inoffsets[ixvel+2], &vz2, SDFtype_sizes[ types[ixvel+2] ]);
        }

        /*dump the outbtab data into the file now-CIE*/
        fwrite(btab, stride, 1, fp);

    }

    /*and we're done! clean up now -CIE*/
    free(members);
    free(addrs);
    free(strides);
    free(nobjs);
    free(starts);
    free(types);
    free(inoffsets);

    free(btab);
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


void renorm(int want[][NISO], float newabund[], float abundarr[], int nparr[], int nnarr[], int Niso, int Nnw, FILE *frp) {
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

        if(is_in_arr == 0) {
            /* figure out where we are */
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

            //printf("index: %3d, isotope: p=%2d n=%2d; put in ",jl,nparr[j],nnarr[j]);

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
    if(fabs(sum - 1.0) > 1.e-3) fprintf(frp, "warning: mass fractions don't sum to 1: %E\n",sum);

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
                   if(newabund[i] < 0.000) fprintf(frp, "warning: negative abundance for %2d %2dafter %d iterations!\n",want[0][i],want[1][i],count);
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
    if(fabs(sum - 1.0) > 1.e-3) fprintf(frp, "renorm warning: mass fractions don't sum to 1: %E\n",sum);
}
