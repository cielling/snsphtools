#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

int main(void) {

    FILE *fp=NULL;
    int i,j;
    float **abundarr, *radbin, *sumabund;
    float totalmass = 0;
    int Nbins, numA;

    /* read in the radial bins - WORKS!!*/
    fp = fopen("inputh.dat", "r");
    if (!fp) printf("error opening inputh.dat\n");

    /*malloc array and read in radial bins. don't know how many, but guess 1000
      and realloc later if we run out of space -CE */
    radbin = (float *)malloc(1000 * sizeof(float));

    Nbins = 0;
    do {
        fscanf(fp, "%21G", &radbin[Nbins++]);
        fscanf(fp, "%*21g"); /* reads in h's. not needed */
        if(!(Nbins % 1000)) realloc(radbin, (Nbins+1000)*sizeof(float));
    } while(radbin[Nbins-1] < 4.3e2);

    fclose(fp);
    fp = NULL;

    /* get Z and N of isotopes - WORKS!!*/
    fp = fopen("abundancereadme", "r");
    if (!fp) printf("error opening file abundancereadme\n");

    /*read in first number, that's the number of isotopes in the file*/
    fscanf(fp, "%3d", &numA);

    fclose(fp);
    fp = NULL;

    /* malloc 2D array for abundances[radial bin][isotope] -CE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    sumabund = (float *)malloc( numA * sizeof(float));
    if (!sumabund) printf("error allocating sumabund\n");
    for( i = 0; i < numA; i++) sumabund[i] = 0.;

    /* read in abundance data - WORKS!!*/
    fp = fopen("abun.dat", "r");
    if (!fp) printf("error opening file abun.dat\n");

    /*read in abundances, one set for each radial bin*/
    for (i=0; i< Nbins; i++){
        totalmass = 0.;
        for (j=0; j < numA; j++) {
            fscanf(fp, "%12G", &abundarr[i][j]);
            sumabund[j] +=abundarr[i][j];
            totalmass += abundarr[i][j];
        }
        printf("sum(row): %.4E\n",totalmass);
    }

    fclose(fp);
    fp = NULL;

    printf("total mass= %.3E grams. H= %.3E  He= %.3E  O= %.3E\n",totalmass, sumabund[numA-1], sumabund[numA-2], sumabund[17]);
    return 0;
}
