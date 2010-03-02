#include<stdlib.h>
#include<stdio.h>

int main(void) {
    FILE *fap;
    int *nparr, *nnarr;
    float *radbin, **abundarr;
    int i, j, numA, Nbins;

    /* now need to read in the abundances, and other necessary stuff */

    /* get Z and N of isotopes - WORKS!!*/
    fap = fopen("abundancereadme", "r");
    if (!fap) printf("error opening file abundancereadme\n");

    /*read in first number, that's the number of isotopes in the file*/
    fscanf(fap, "%3d", &numA);

    nparr = (int *)malloc(numA * sizeof(int));
    nnarr = (int *)malloc(numA * sizeof(int));

    for ( i = 0; i < numA; i++) fscanf(fap, "%d", &nparr[i]);
    for ( i = 0; i < numA; i++) fscanf(fap, "%d", &nnarr[i]);

    fclose(fap);
    fap = NULL;


    /* read in the radial bins - WORKS!!*/
    fap = fopen("inputh.dat", "r");
    if (!fap) printf("error opening inputh.dat\n");

    /*malloc array and read in radial bins. don't know how many, but guess 1000
      and realloc later if we run out of space -CE */
    radbin = (float *)malloc(1000 * sizeof(float));

    Nbins = 0;
    do {
        fscanf(fap, "%21G", &radbin[Nbins++]);
        fscanf(fap, "%*21g"); /* reads in h's. not needed */
        if(!(Nbins % 1000)) realloc(radbin, (Nbins+1000)*sizeof(float));
    } while(radbin[Nbins-1] < 4.3e2);

    fclose(fap);
    fap = NULL;

    /* malloc 2D array for abundances[radial bin][isotope] -CE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    /* read in abundance data - WORKS!!*/
    fap = fopen("abun.dat", "r");
    if (!fap) printf("error opening file abun.dat\n");

    /*read in abundances, one set for each radial bin*/
    for (i=0; i< Nbins; i++){
        for (j=0; j < numA; j++) {
            fscanf(fap, "%12G", &abundarr[i][j]);
        }
    }

    fap = fopen("abund.txt", "w");
    if (!fap) printf("error opening abund.txt\n");
    for( i = 0; i < Nbins; i++ )
        fprintf(fap, "%.6G  %.6G\n", radbin[i], abundarr[i][0]);

    fclose(fap);

    return 0;
}
