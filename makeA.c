#include<stdlib.h>
#include<stdio.h>
#include<string.h>

int main(int argc, char *argv[]) {
    FILE *fp;
    int *nparr, *nnarr;
    float *radbin, **abundarr, sumX = 0.0;
    char **isoarr, myline[2000];
    int i, j, numA, Nbins;
    int calc_he;

    if(argc != 2)
       printf("Error, need value for calc_He_flag!\n");

    calc_he = atoi(argv[1]);

    if(calc_he) printf("calculating He abundance\n");

    /* now need to read in the abundances, and other necessary stuff */

    /* get Z and N of isotopes - WORKS!!*/
    fp = fopen("abundances-README", "r");
    if (!fp) printf("error opening file abundancereadme\n");

    /*read in first number, that's the number of isotopes in the file*/
    //fscanf(fp, "%d", &numA);
    numA = 177;

    nparr = (int *)malloc(numA * sizeof(int));
    nnarr = (int *)malloc(numA * sizeof(int));
    isoarr = (char **)malloc(numA * sizeof(char *));
    for( i = 0; i < numA; i++)
        isoarr[i] = (char *)malloc( 5*sizeof(char));

/*
    for ( i = 0; i < numA; i++) fscanf(fp, "%d", &nparr[i]);
    for ( i = 0; i < numA; i++) fscanf(fp, "%d", &nnarr[i]);
    for ( i = 0; i < numA; i++) fscanf(fp, "%s", &isoarr[i]);
*/
    for( i = 0; i < numA; i++) {
        fscanf(fp, "%*d %s %d %d", isoarr[i], &nparr[i], &nnarr[i]);
    }

    fclose(fp);
    fp = NULL;


    /* read in the radial bins - WORKS!!*/
    fp = fopen("inputh.dat", "r");
    if (!fp) printf("error opening inputh.dat\n");

    do {
        fgets(myline, 2000, fp);
        Nbins++;
    } while(!feof(fp));

    fclose(fp);
    fp = NULL;
    fp = fopen("inputh.dat", "r");
    /*malloc array and read in radial bins. don't know how many, but guess 1000
      and realloc later if we run out of space -CE */
    radbin = (float *)malloc(Nbins * sizeof(float));

    Nbins = 0;
    do {
        fscanf(fp, "%21G", &radbin[Nbins++]);
        fscanf(fp, "%*21g"); /* reads in h's. not needed */
    } while(!feof(fp));

    fclose(fp);
    fp = NULL;

    Nbins--;

    /* malloc 2D array for abundances[radial bin][isotope] -CE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    /* read in abundance data - WORKS!!*/
    fp = fopen("abundances", "r");
    if (!fp) printf("error opening file abun.dat\n");

    /*read in abundances, one set for each radial bin*/
    for (i=0; i< Nbins; i++){
        sumX = 0.0;
        for (j=0; j < numA; j++) {
            fscanf(fp, "%12G", &abundarr[i][j]);
            if (calc_he)
                sumX += abundarr[i][j]; /* add ALL X_i, incl. He */
        }
        if(calc_he) {
            sumX -= abundarr[i][numA-1]; /* subtract wrong X(He) */
            abundarr[i][numA-1] = ( (1.0 -sumX) < 0.0 ? 0.0 : (1.0-sumX) );
        }
    }

    fp = fopen("abund.txt", "w");
    if (!fp) printf("error opening abund.txt\n");

    fprintf(fp, "%d\n%d\n", numA, Nbins);

    for ( i = 0; i < numA; i++ )
        fprintf(fp, "%d ", nparr[i]);
    fprintf(fp, "\n");
    for ( i = 0; i < numA; i++ )
        fprintf(fp, "%d ", nnarr[i]);
    fprintf(fp, "\nr ");
    for ( i = 0; i < numA; i++ )
        fprintf(fp, "%s ", isoarr[i]);
    fprintf(fp, "\n");

    for( i = 0; i < Nbins; i++ ) {
        fprintf(fp, "%12.6E ", radbin[i]);
        for( j = 0; j < numA; j++){
    /*        if(nparr[j] == nnarr[j])*/
            fprintf(fp, "%12.6E ", abundarr[i][j]);
        }
        fprintf(fp, "\n");
    }
/*
        fprintf(fp, "%.6G  %.6G\n", radbin[i], abundarr[i][0]);
*/

    fclose(fp);

    return 0;
}
