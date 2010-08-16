/*
   PURPOSE:
	to add an abundance table to an initial conditions sdf file.
        This version should also work for large files, since it reads in the
        sdf file line-by-line.

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

static const int NISO = 22;

int main(int argc, char *argv[])
{
    FILE *fp = NULL;
    FILE *fap = NULL;
    int i, j, k, nvecs, ju, jm, jl;
    int numA, Nbins, len, counter, flag = 0;
    int *nparr, *nnarr;
    int want[2][NISO]={{0,1,2,6,7,8,10,12,14,15,16,18,20,20,22,24,26,26,26,27,28,28},/*Z*/
                       {1,0,2,6,7,8,10,12,14,16,16,18,20,24,22,24,26,30,32,29,28,30}}; /*A-Z*/
    char elements[33][3] = {{"N"},{"h"},{"he"},{"li"},{"be"},{"b"},{"c"},{"n"},{"o"},{"f"},{"ne"},
                            {"na"},{"mg"},{"al"},{"si"},{"p"},{"s"},{"cl"},{"ar"},{"k"},{"ca"},
                            {"sc"},{"ti"},{"v"},{"cr"},{"mn"},{"fe"},{"co"},{"ni"},{"cu"},
                            {"zn"},{"ga"},{"ge"}};
    char mychar[5];
    size_t stride = 0, outstride = 0;
    float **abundarr; /*should this be void * ? */
    double x, y, z, radius;
    float *radbin;


/*malloc memory space for the respective features of the struct-CE*/
    //for(j=0;j<32;j++) elements[j] = (char *)malloc(3 * sizeof(char));

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
            fscanf(fap, "%12E", &abundarr[i][j]);
        }
    }

    fclose(fap);
    fap = NULL;

    fp = fopen("run3g_50_abun.dat", "w");
    if(!fp) printf("error opening file 'rung3g_50_abun.dat'\n");

    fprintf(fp,"%d\n",numA);
    fprintf(fp,"%d\n",Nbins);
    for (j=0; j < numA; j++) fprintf(fp,"%d ",nparr[j]);
    fprintf(fp,"\n");
    for (j=0; j < numA; j++) fprintf(fp,"%d ",nnarr[j]);
    fprintf(fp,"\n");

    fprintf(fp,"%s ","r");
    for (j=0; j < numA; j++) fprintf(fp,"%s%-d ",elements[ nparr[j] ], nparr[j]+nnarr[j]);
    fprintf(fp,"\n");

    for(i=0;i<Nbins;i++) {
        fprintf(fp,"%10E ",radbin[i]);
        for(j=0;j<numA;j++) {
            fprintf(fp,"%10E ", (double)abundarr[i][j]);
        }
        fprintf(fp,"\n");
    }

    return 0;
}
