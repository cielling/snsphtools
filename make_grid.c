#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<errno.h>
#include<SDF.h>
#include"physics_sph.h"

//const int NDIM=3;

void initargs(int argc, char* argv[], FILE** fp, int *nx, int *ny, int *nz);

int main(int argc, char* argv[])
{
	FILE *fp=NULL;
	SPHoutbody *outbtab, OneBody;
	int i, j, k, index;
	int nx, ny, nz;
	float nx2, ny2, nz2;
	int ngrid;
	float dx, dy, dz, dh;
	float ii, jj, kk;
	float R0=4.e5;

	/*initargs: check if file exists, ask if overwrite*/
	initargs(argc, argv, &fp, &nx, &ny, &nz);

	printf("set nx= %d ny= %d nz= %d\n",nx, ny, nz);

	/*create the grid */
	ngrid = nx*ny*nz;
	outbtab = malloc( ngrid*sizeof(SPHoutbody) );
	if(!outbtab) fprintf(stderr, "error allocating outbtab!\n");
	
	dx = 2.*R0/(float)(nx-1);
	dy = 2.*R0/(float)(ny-1);
	dz = 2.*R0/(float)(nz-1);
	nx2 = (float)(nx-1)*0.5;
	ny2 = (float)(ny-1)*0.5;
	nz2 = (float)(nz-1)*0.5;

	dh = (dx+dy+dz)/3. *0.5;

	/* try different binning also, like log, sqrt */

	for( i=0; i<nx; ++i) {

		ii = (float)i-nx2;
		for( j=0; j<ny; ++j) {

			jj = (float)j-ny2;
			for( k=0; k<nz; ++k) {

				kk = (float)k-nz2;
				index = i*ny*nz + j*nz + k;
				
				outbtab[index].pos[0] = (double)ii*dx;
				outbtab[index].pos[1] = (double)jj*dy;
				outbtab[index].pos[2] = (double)kk*dz;
				outbtab[index].mass = 1.;
				outbtab[index].h = dh;
				outbtab[index].rho = 0.;
				/*
				outbtab[index].nterms = 0.;
				outbtab[index].temp = 0.;
				outbtab[index].vel[0] = 0.;
				outbtab[index].vel[1] = 0.;
				outbtab[index].vel[2] = 0.;

				OneBody.pos[0] = (double)ii*dx;
				OneBody.pos[1] = (double)jj*dy;
				OneBody.pos[2] = (double)kk*dz;
				OneBody.mass = 0.;
				OneBody.h = dh;

				outbtab[index].pos[0] = OneBody.pos[0];
				outbtab[index].pos[1] = OneBody.pos[1];
				outbtab[index].pos[2] = OneBody.pos[2];
				outbtab[index].mass = OneBody.mass;
				outbtab[index].h = OneBody.h;
				*/

			}
		}
	}
	
	/*write the SDF header */
	/*
	fprintf(fp, "# SDF\n");
	fprintf(fp, "parameter byteorder=%#x;\n", SDFcpubyteorder() );
	fprintf(fp, "int npart = %d;\n", ngrid);
	fprintf(fp, "int ndim = %d;\n", 3);
	fprintf(fp, "float R0 = %f;\n", R0);
	fprintf(fp, "struct {\n");
	fprintf(fp, "\tdouble x, y, z;\n");
	fprintf(fp, "\tfloat h;\n");
	fprintf(fp, "\tfloat mass;\n");
	fprintf(fp, "\tfloat rho;\n");
	fprintf(fp, "}[%d]\n",ngrid);
	fprintf(fp, "\#\^L\n\# SDF-EOH\n");
	*/
		/*
		*/
	SDFwrite(argv[1],ngrid, ngrid, outbtab, sizeof(SPHoutbody), SPHOUTBODYDESC,
		"npart", SDF_INT, ngrid,
		"iter", SDF_INT, 0,
		"dt", SDF_FLOAT, 0.,
		"eps", SDF_FLOAT, 0.,
		"Gnewt", SDF_FLOAT, 0.,
		"tolerance", SDF_FLOAT, 0.,
		"frac_tolerance", SDF_FLOAT, 0.,
		"ndim", SDF_INT, 3,
		"tpos", SDF_FLOAT, 0.,
		"tvel", SDF_FLOAT, 0.,
		"gamma", SDF_FLOAT, 1.67,
		"R0", SDF_FLOAT, R0,
		NULL);
		

	/*write the SDF body */

	fclose(fp);
}


void initargs(int argc, char* argv[], FILE** fp, int *nx, int *ny, int *nz) {

	char input;

	if( argc != 5 ) {
		fprintf(stderr, "Usage: %s outfile nx ny nz\n", argv[0]); 
		exit(1);
	}


	if( access(argv[1], F_OK) == 0 ){

		fprintf(stderr, "%s: %s exists, overwrite? y/n\n", argv[0], argv[1]);
		input = getc(stdin);

		if((input != 'y') && (input != 'Y')) {
			fprintf(stderr, "exiting.\n");
			exit(1);
		}
	}

	*fp = fopen(argv[1], "w");
	
	if( fp == NULL ) {
		fprintf(stderr, "Error, could not open file!\n");
		fprintf(stderr, "%s: %s\n", argv[0], strerror(errno));
		exit(errno);
	}

	*nx = atoi(argv[2]);
	*ny = atoi(argv[3]);
	*nz = atoi(argv[4]);

}
