/*
   PURPOSE:
	to get the yield of a list of isotopes

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
#include <math.h>
#include <SDF.h>

typedef enum SDF_type_enum SDF_type_t;

typedef union {
    int i;
    float f;
    double d;
} datum_t;

static int **isotope;
static int countis;

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp, FILE **fp2);
static void writestructs(SDF *sdfp, FILE *fp, FILE *fp2, int *partids, int npartids, float *ppyields);

int	read_se_yields(char *argv[], int **partids, float **ppyields, FILE *fp2);
int make_spec_names(char ***specarr, char spec, int num);

int main(int argc, char *argv[])
{
    SDF *sdfp = NULL;
    FILE *fp = NULL, *fp2 = NULL;
	int *partids, npartids;
	float *ppyields;

    initargs(argc, argv, &sdfp, &fp, &fp2);

	npartids = read_se_yields(argv, &partids, &ppyields, fp2);
	if( npartids <= 0)
		printf("WARNING! no post-proc'd particles found!\n");

	fclose(fp2);
	fp2 = fopen(argv[3], "r");

    writestructs(sdfp, fp, fp2, partids, npartids, ppyields);

    fclose(fp);
    fclose(fp2);
    SDFclose(sdfp);

	free(partids);

    return 0;
}

static void initargs(int argc, char *argv[], SDF **sdfp, FILE **fp, FILE **fp2)
{
    char input;

    if (argc != 5) {
	fprintf(stderr, "Usage: %s SDFfile outfile ctl-file postproc-file\n", argv[0]);
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

    *fp2 = fopen(argv[3], "r");
    if (*fp2 == NULL) {
	fprintf(stderr, "%s: %s\n", argv[3], strerror(errno));
	exit(errno);
    }

}

int read_se_yields(char *argv[], int **partid, float **ppyields, FILE *fp2)
{
	FILE *fp=NULL;
	int i, j, nlines=0;
	char temp[21], line[101];
	int counti=0, nn, np;
	float mass;


	/* read get.isotopes-ctl file */
	while( !feof(fp2) ){
		fgets(line, 100, fp2);
		counti++;
	}
	rewind(fp2);
	counti--;

	isotope = (int **)malloc( 2*sizeof(int *) );
	isotope[0] = (int *)malloc( counti*sizeof(int) );
	isotope[1] = (int *)malloc( counti*sizeof(int) );

	*ppyields = (float *)malloc( counti*sizeof(float) );


	counti = 0;

	while( !feof(fp2) ){
		fscanf(fp2, "%d %d", &isotope[0][counti], &isotope[1][counti]);
		counti++;
	}
	counti--;
	countis = counti;

	rewind(fp2);


	/* read post-processed particles file */
	fp = fopen(argv[4], "r");
	if(fp == NULL) {
		printf("error opening file! %s\n",argv[4]);
		return 0;
	}


	while( !feof(fp) ) {
		fgets(line, 100, fp);
		nlines++;
	}

	rewind(fp);

	*partid = (int *)malloc( nlines*sizeof(int) );
	fgets(line, 100, fp);

	for(i = 0; i < nlines; i++) {
		fscanf(fp, "%s %d %*s %*s %*e", temp, &((*partid)[i]));
		if( strncmp(temp, "nn", strlen(temp)) == 0)
			break;
<<<<<<< HEAD
	}
=======
        }
>>>>>>> snsphtools-share/master
	fgets(line, 100, fp);

	/* for now, just assume that we never want to get the yield of n=1,p=1 */
	while( !feof(fp) ) {
		fscanf(fp, "%*s%*s %d %*s %*s %d %*s %*s %e %*s", &nn, &np,&mass);
		for( j = 0; j < counti; j++) {
			if((nn == isotope[1][j]) && (np == isotope[0][j])) {
			(*ppyields)[j] = mass/1.99e33;
			}
		}
	}

/*
	free(isotope[1]);
	free(isotope[0]);
	free(isotope);
*/

	fclose(fp);

	return i;
}


/*this writes the actual data*/
static void writestructs(SDF *sdfp, FILE *fp, FILE *fp2, int *partids, int npartids, float *ppyields)
{
    int i, j, k, nvecs, nmembers;
	int ii, postprocd;
    char **vecs, **members;
	char **pnames, **nnames;
	char line[101];
    SDF_type_t *types, type;
    size_t stride = 0;
    void *btab;
    void **addrs;
    int *inoffsets, *lines, *strides, *starts;
    int INCR=1, flag=0;
    int nlines = 1, nrecs;
    int index[6], counti, ident, identmin;
	int n_iso, **iso_arr, *gotcha;
    float rho, *yield, mass, *Xel, massCF, temp, total;
	int **isotop;
    /*make INCR and nlines user input */

/* this does not attempt to load the whole file into memory, does it? -CE */
    nvecs = SDFnvecs(sdfp);
    vecs = SDFvecnames(sdfp);

    SDFgetint(sdfp, "npart", &nrecs);

	/* prompt user for starting particle id */
	printf("enter starting particle id; '-1' for all particles\n");
	scanf("%d", &identmin);

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
        if (strncmp(vecs[i], "mass", strlen(vecs[i])) == 0) index[1]=i;
        if (strncmp(vecs[i], "f1", strlen(vecs[i])) == 0) index[2]=i;
        if (strncmp(vecs[i], "p1", strlen(vecs[i])) == 0) index[3]=i;
        if (strncmp(vecs[i], "ident", strlen(vecs[i])) == 0) index[4]=i;
        if (strncmp(vecs[i], "rho", strlen(vecs[i])) == 0) index[5]=i;

		if (flag) ++nmembers;
    }
    printf("nmembers = %d\n",nmembers);

<<<<<<< HEAD
	n_iso = index[3]-index[2];
	//n_iso = 20;
=======
	//n_iso = index[3]-index[2];
	n_iso = 20;
>>>>>>> snsphtools-share/master
	printf("%d\n", n_iso);
	iso_arr = (int**)malloc( 2*sizeof(int *) );
	iso_arr[0] = (int *)malloc( n_iso*sizeof(int) );
	iso_arr[1] = (int *)malloc( n_iso*sizeof(int) );

<<<<<<< HEAD
    make_spec_names(&nnames, 'm', n_iso);
=======
    make_spec_names(&nnames, 'n', n_iso);
>>>>>>> snsphtools-share/master
    make_spec_names(&pnames, 'p', n_iso);
	printf("WARNING: make sure you've got 'm' vs 'n' and n_iso correct!\n");

	if( SDFgetfloat(sdfp, "massCF", &massCF) != 0) massCF = 1.998e27;
	massCF = massCF/1.998e33;

	for( i=0; i < n_iso; i++) {
		SDFgetint(sdfp,pnames[i], &iso_arr[0][i]);
		SDFgetint(sdfp,nnames[i], &iso_arr[1][i]);
		printf("%3s = %2d\t", pnames[i], iso_arr[0][i]);
		printf("%3s = %2d\n", nnames[i], iso_arr[1][i]);
	}

	counti = 0;
	while( !feof(fp2) ){
		fgets(line, 100, fp2);
		counti++;
	}

	printf("getting %d isotopes\n",--counti);
	//counti = countis+1;

	isotop = (int **)malloc( 2*sizeof(int *) );
	isotop[0] = (int *)malloc( counti*sizeof(int) );
	isotop[1] = (int *)malloc( counti*sizeof(int) );

	rewind(fp2);

	gotcha = (int *)malloc((counti+2)*sizeof(int));
	Xel = (float *)malloc( counti*sizeof(float) );
	yield = (float *)malloc( counti*sizeof(float) );

	counti=0;
	ii = 0;

	while( !feof(fp2) ){
		fscanf(fp2, "%d %d", &isotop[0][ii], &isotop[1][ii]);
		//fscanf(fp2, "%*d %*d");

		for(k = 0; k < n_iso; k++) {
			/* find isotope */
<<<<<<< HEAD
				if( (iso_arr[0][k] == isotop[0][ii]) &&
=======
			if( (iso_arr[0][k] == isotop[0][ii]) &&
>>>>>>> snsphtools-share/master
				(iso_arr[1][k] == isotop[1][ii]) ) {
				gotcha[counti++] = k + index[2];
				break;
			}
		}
		ii++;
	}

	//--counti;
	printf("counti: %d\n",counti);

	/* sort */
	if(counti>1) {

	for(i = 0; i < counti; i++) {
		for(j = 0; j < counti-1; j++ ) {
			if(gotcha[j] > gotcha[j+1]) {
				temp = gotcha[j+1];
				gotcha[j+1] = gotcha[j];
				gotcha[j] = temp;
			}
		}
	}

	}

	for( i = 0; i < counti; i++) {
		printf("%d\t",gotcha[i]);
	}
	printf("\n");

/*malloc memory space for the respective features of the struct-CE*/
    members = (char **)malloc((counti+3) * sizeof(char *));
    addrs = (void **)malloc((counti+3)* sizeof(void *));
    types = (SDF_type_t *)malloc((counti+3) * sizeof(SDF_type_t));
    inoffsets = (int *)malloc((counti+3)* sizeof(int));
    strides = (int *)malloc((counti+3) * sizeof(int));
    lines = (int * )malloc((counti+3) * sizeof(int));
    starts = (int * )malloc((counti+3) * sizeof(int));

/*one by one, go through the fields in the column, i.e. members of the struct?-CE*/

    flag = 0;
    for (i = 0, stride = 0, j = 0; i < nvecs; ++i) {
<<<<<<< HEAD
		if(i == index[1]) {
=======
            if(i == index[1]) {
>>>>>>> snsphtools-share/master
			/* get mass */
			members[0] = vecs[i];
            lines[0] = 1;//INCR;
            starts[0] = 0;
            types[0] = SDFtype(members[0], sdfp);
            inoffsets[0] = stride;
            stride += SDFtype_sizes[ SDFtype(vecs[i],sdfp) ];/*applies to output*/
            printf("member = %s offset: %d \n", members[0],inoffsets[0]);
		}
		if(i == gotcha[j] && j<counti) {
            members[j+1] = vecs[i];
            lines[j+1] = 1;//INCR;
            starts[j+1] = 0;
            types[j+1] = SDFtype(members[j+1], sdfp);
            inoffsets[j+1] = stride;
            stride += SDFtype_sizes[ SDFtype(vecs[i],sdfp) ];
            printf("member = %s offset: %d \n",members[j+1],inoffsets[j+1]);
			j++;
		}
    }

	/* get particle id */
    members[counti+1] = vecs[index[4]];
    lines[counti+1] = 1;//INCR;
    starts[counti+1] = 0;
    types[counti+1] = SDFtype(members[counti+1], sdfp);
    inoffsets[counti+1] = stride;
    stride += SDFtype_sizes[ SDFtype(vecs[index[4]],sdfp) ];
    printf("member = %s offset: %d \n",members[counti+1],inoffsets[counti+1]);

	/* get particle density */
    members[counti+2] = vecs[index[5]];
    lines[counti+2] = 1;//INCR;
    starts[counti+2] = 0;
    types[counti+2] = SDFtype(members[counti+2], sdfp);
    inoffsets[counti+2] = stride;
    stride += SDFtype_sizes[ SDFtype(vecs[index[5]],sdfp) ];
    printf("member = %s offset: %d \n",members[counti+2],inoffsets[counti+2]);

	nmembers = j+3;

    btab = (void *)malloc(stride * nlines);

    /*calculate the byte offset in memory to the address of the next member-CE*/
	addrs[0] = (char *)btab;
    strides[0] = stride;		/* how did this run correctly without this line??? */
	for ( i = 1; i < nmembers; i++) {
            addrs[i] = addrs[i-1] + SDFtype_sizes[ types[i-1] ];
            strides[i] = stride;
        }


    printf("reading in %d lines ...\n",nrecs);


	ii = 0;

<<<<<<< HEAD
		/*try reading in one line at a time */
		for( j = 0; j < nrecs; j++) {
			for( i = 0; i < nmembers; i++)
				starts[i] = j;

			if(!(j%50000)) printf("iter: %d\n",j);
			SDFseekrdvecsarr(sdfp, nmembers, members, starts, lines, addrs, strides);
			/*temporary solution*/

			mass = *((float *)(btab + inoffsets[0]));
			ident = *((int *)(btab + inoffsets[counti+1]));
			rho = *((float *)(btab + inoffsets[counti+2]));
			if(!(j%50000)) printf("rho: %e\n",rho);
=======
        /*try reading in one line at a time */
        for( j = 0; j < nrecs; j++) {
                for( i = 0; i < nmembers; i++)
                starts[i] = j;

		if(!(j%50000)) printf("iter: %d\n",j);
                SDFseekrdvecsarr(sdfp, nmembers, members, starts, lines, addrs, strides);
			/*temporary solution*/

                mass = *((float *)(btab + inoffsets[0]));
                ident = *((int *)(btab + inoffsets[counti+1]));
                rho = *((float *)(btab + inoffsets[counti+2]));
                if(!(j%50000)) printf("rho: %e\n",rho);
>>>>>>> snsphtools-share/master

			postprocd = 0; /* only add abundances of non-post proc'd particles */

			if( npartids > 0 ) {
			if( (partids[ii] == ident) && (ii < npartids)) {
				postprocd = 1;
				ii++;
			} else {
				postprocd = 0;
			}
			}

			if( postprocd==0 ) {
				for(k = 0; k< counti; k++) {
<<<<<<< HEAD
					Xel[k] = *((float *)(btab + (k+1)*sizeof(float) ));
						if( rho >= 5.0e-3) {
						yield[k] += Xel[k] * mass * massCF;
                    }
=======
                                Xel[k] = *((float *)(btab + (k+1)*sizeof(float) ));
                                if( rho >= 5.0e-3) {
					yield[k] += Xel[k] * mass * massCF;
                                }
>>>>>>> snsphtools-share/master
				}
			}

		}

		printf("got mass?\n");
		total = 0;

		fprintf(fp, "%4s %4s %14s\t%14s\t%14s\n", "p", "n", "un-proc'd", "proc'd", "total");
		for(k = 0; k < counti; k++) {
			fprintf(fp, "%4d %4d %14.6e", iso_arr[0][ gotcha[k]-index[2] ],
					iso_arr[1][ gotcha[k]-index[2] ], yield[k]);
			if(!(yield[k]!=yield[k]))
				total += yield[k];
			for(ii = 0; ii < counti; ii++) {
				if((iso_arr[0][ gotcha[k]-index[2] ] == isotop[0][ii]) &&
					(iso_arr[1][ gotcha[k]-index[2] ] == isotop[1][ii]) ) {
					fprintf(fp, "\t%14.6e\t%14.6e\n", ppyields[ii], (yield[k]+ppyields[ii]));
				}
			}
		}
		fprintf(fp, "\n\nunprocessed mass: %e\n",total);

		printf("total mass: %e\n",total);

/*and we're done! clean up now -CE: if it ever works*/
/*    free(members);
    free(btab);
    free(addrs);
*/
    /*free(outbtab);*/
	free(lines);
    free(types);
    free(inoffsets);
	free(strides);
	free(iso_arr[0]);
	free(iso_arr[1]);
	free(iso_arr);
	free(yield);
	free(Xel);
	free(gotcha);

}






int make_spec_names(char ***specarr, char spec, int num)
{
    int i;
    char tmpchr[20];

    *specarr = (char **)malloc(num * sizeof(char *) );

    for( i = 0; i < num; i++ ){

        sprintf( tmpchr, "%c%-d\0", spec, (i+1));
        *(*specarr + i) = (char *)malloc( ( strlen(tmpchr) + 1) * sizeof(char) );
        sprintf((*specarr)[i], "%s",tmpchr);
        //printf("isotope specifier: %s  %d\n", specarr[i],(int)strlen(specarr[i]));

    }
    return i;
}
