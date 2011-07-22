/*

   PURPOSE:
	to read the bndry information from the header of a sequence of
	sdf dump files and compile into a table. Can read either 4- or
	5- digit long extensions of the files (i.e. the .<step> extension)
	but not both at the same time.

   NOTE:
	apparently, the SDF* functions can't read sdf files for which one
	or more of the header quantities equal 'nan'. ....
	This is so far the reason for the 'error opening sdf file'
	messages.

   INPUT FILES:
	<input-ctl-file>
		this files drives the program. It contains:
		- the sdf file base name, including absolute path and the '.'
		  before the step number extension
		- number of first file in sequence
		- number of last file in sequence
		- the step size or increment between file numbers
		- name of the output file

   COMPILE:
	with Makefile. un-comment appropriate lines.
	This routine needs some libraries from the tree code and SDF routines,
	so make sure that TREEHOME is set, and the tree code (SNSPH) compiled
	once for serial use (i.e. without PAROS flag).

   SYNTAX:
	get_bndry <input-ctl-file>

   METHOD:
	@ open <input-clt-file> to determine sequence of sdf files
	@ for each sdf file, read desired quantities from header
	@ write the first line of the read in data block to designated output file

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

typedef struct {
    float pos[3];
    float h;
    float mass;
    float rho;
} SPHbody;

//static float **bndrydata;

int INCR = 25;

static void writeinit(char *outfile, FILE **fp, int first, int last, int dn, float ***bndrydata);
static void writescalars(SDF *sdfp, FILE *fp);
static void writestructs(SDF *sdfp, FILE *fp);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL;
	FILE *fctl = NULL;
	float test;
    int n, dgtl;
	int first, last, dn;
	char sdfbase[200], outfile[50], *sdfname;
	float **bndrydata;

	/* read in prep stuff */
	fctl = fopen(argv[1], "r");
	if(!fctl) printf("error opening ctl-bndry.dat\n");

	fscanf(fctl, "%s", sdfbase);
	fscanf(fctl, "%d", &first);
	fscanf(fctl, "%d", &last);
	fscanf(fctl, "%d", &dn);
	fscanf(fctl, "%s", outfile);

	fclose(fctl);

	/* case of step extension being 4 vs 5 digits long.
	   NOTE: only handles either or at this point! ~CIE
	*/
	if( last >= 10000)
		dgtl = 6;
	else
		dgtl = 5;


	sdfname = (char *)malloc( (strlen(sdfbase)+dgtl)*sizeof(char) );
	strcpy(sdfname, sdfbase);

    writeinit(outfile, &fp, first, last, dn, &bndrydata);

	for( n=first; n<=last; n+=dn) {

	/* create filename, open file */
	if(dgtl == 5)
		sprintf(&sdfname[ strlen(sdfbase) ], "%04d", n);
	else
		sprintf(&sdfname[ strlen(sdfbase) ], "%05d", n);

	printf("%s\n", sdfname);

	//sdfp = SDFopen(NULL, "/scratch/cellinge/runsnsph/r3g_1M_cco_sph.2000");
	sdfp = SDFopen(NULL, sdfname);
	if(!sdfp) printf("error opening sdf file %s!\n",sdfname);

    writescalars(sdfp, fp);/*writes the header for the scalars (non-structs)*/

    /*fclose(fp);*/
    SDFclose(sdfp);

	}

	fclose(fp);

    return 0;
}


static void writeinit(char *outfile, FILE **fp, int first, int last, int dn, float ***bndrydata)
{
	int i, j, num;


    *fp = fopen(outfile, "w");
    if (*fp == NULL) {
	fprintf(stderr, "%s: %s\n", outfile, strerror(errno));
	exit(errno);
    }

	num = (int)((last-first)/dn);
	*bndrydata = (float **)malloc( num*sizeof(float*) );
	if(*bndrydata == NULL) printf("allocation error: bndrydata\n");
	for( i=0; i<num; i++) {
		(*bndrydata)[i] = (float *)malloc( 11*sizeof(float) );
		if((*bndrydata)[i] == NULL) printf("allocation error:bndrydata[i]\n");
		for( j=0; j<11; j++)
			(*bndrydata)[i][j] = -2.0;
	}

	fprintf(*fp, "%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
			"t", "x", "y", "z", "px", "py", "pz", "lx", "ly", "lz", "mass");

}

static void writescalars(SDF *sdfp, FILE *fp)
{
    int i, nvecs, nrecs;
    int num=0;
    SDF_type_t type;
	float bndrydata[1][11];

    nvecs = SDFnvecs(sdfp);/*figure out number of lines in the header, basically-CE*/

		/* for some reason, using the bndrydata that is allocated in
		   writeinit doesn't work here; I'm getting seg faults in the
		   fprintf statement from it below, the SDFget* work with it
		   though ....
		*/
		SDFgetfloat(sdfp, "tpos", &bndrydata[num][0]);
		SDFgetfloat(sdfp, "bndry_x", &bndrydata[num][1]);
		SDFgetfloat(sdfp, "bndry_y", &bndrydata[num][2]);
		SDFgetfloat(sdfp, "bndry_z", &bndrydata[num][3]);
		SDFgetfloat(sdfp, "bndry_px", &bndrydata[num][4]);
		SDFgetfloat(sdfp, "bndry_py", &bndrydata[num][5]);
		SDFgetfloat(sdfp, "bndry_pz", &bndrydata[num][6]);
		SDFgetfloat(sdfp, "bndry_lx", &bndrydata[num][7]);
		SDFgetfloat(sdfp, "bndry_ly", &bndrydata[num][8]);
		SDFgetfloat(sdfp, "bndry_lz", &bndrydata[num][9]);
		SDFgetfloat(sdfp, "bndry_mass", &bndrydata[num][10]);

	for( i=0; i<11; i++)
		fprintf(fp, "%12.4e", bndrydata[num][i]);
	fprintf(fp, "\n");



}

