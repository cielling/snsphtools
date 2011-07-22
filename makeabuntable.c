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
    FILE *afp = NULL;
    int i, j, k, nvecs, ju, jm, jl, m, n;
    int numA, Nbins, Abins = 0, NNW, len, counter, flag = 0;
    int *nparr, *nnarr, **nwlist;
    int want[2][NISO]={{0,1,2,6,7,8,10,12,14,15,16,18,20,20,22,24,26,26,26,27,28,28},/*Z*/
                       {1,0,2,6,7,8,10,12,14,16,16,18,20,24,22,24,26,30,32,29,28,30}}; /*A-Z*/
    char elements[33][3] = {{"N"},{"h"},{"he"},{"li"},{"be"},{"b"},{"c"},{"n"},{"o"},{"f"},{"ne"},
                            {"na"},{"mg"},{"al"},{"si"},{"p"},{"s"},{"cl"},{"ar"},{"k"},{"ca"},
                            {"sc"},{"ti"},{"v"},{"cr"},{"mn"},{"fe"},{"co"},{"ni"},{"cu"},
                            {"zn"},{"ga"},{"ge"}};
    char mychar[5],tmpchr[40];
    size_t stride = 0, outstride = 0;
    float **abundarr; /*should this be void * ? */
    double *He;
    double x, y, z, radius;
    float *radbin;


/*malloc memory space for the respective features of the struct-CE*/
    //for(j=0;j<32;j++) elements[j] = (char *)malloc(3 * sizeof(char));

    /* now need to read in the abundances, and other necessary stuff */

    /* get Z and N of isotopes - WORKS!!*/
    afp = fopen("abundancereadme", "r");
    if (!afp) printf("error opening file abundancereadme\n");

    /*read in first number, that's the number of isotopes in the file*/
    fscanf(afp, "%3d", &numA);

    nparr = (int *)malloc(numA * sizeof(int));
    nnarr = (int *)malloc(numA * sizeof(int));

    for ( i = 0; i < numA; i++) fscanf(afp, "%d", &nparr[i]);
    for ( i = 0; i < numA; i++) fscanf(afp, "%d", &nnarr[i]);

    fclose(afp);
    afp = NULL;


    /* read in the radial bins - WORKS!!*/
    afp = fopen("inputh.dat", "r");
    if (!afp) printf("error opening inputh.dat\n");

    /*malloc array and read in radial bins. don't know how many, but guess 1000
      and realloc later if we run out of space -CE */
    radbin = (float *)malloc(1000 * sizeof(float));

    Nbins = 0;
    do {
        fscanf(afp, "%21G", &radbin[Nbins++]);
        fscanf(afp, "%*21g"); /* reads in h's. not needed */
        if(!(Nbins % 1000)) realloc(radbin, (Nbins+1000)*sizeof(float));
    } while(!feof(afp));
    /*} while(radbin[Nbins-1] < 4.3e2);*/

    fclose(afp);
    afp = NULL;

    /* malloc 2D array for abundances[radial bin][isotope] -CE */
    abundarr = (float **)malloc(Nbins * sizeof( float *));
    if (!abundarr) printf("error allocating abundarr\n");

    for (i=0; i < Nbins; i++) {
        abundarr[i] = (float *)malloc( numA * sizeof(float));
        if (!abundarr[i]) printf("error allocating abundarr[i]\n");
    }

    /* read in abundance data - WORKS!!*/
    afp = fopen("abundances", "r");
    if (!afp) printf("error opening file abun.dat\n");

    He = (double *)malloc(Nbins * sizeof(double));

    /*read in abundances, one set for each radial bin*/
    i=0;
    do {
        He[i] = 0.0;
        for(j=0; j<numA; j++) {
            fscanf(afp, "%12E", &abundarr[i][j]);
            He[i] += (double)abundarr[i][j];
        }
        (He[i] > 1.0) ? (abundarr[i][numA-1] = 0.0) :
            (abundarr[i][numA-1] = (float)((double)1.0 - He[i]));
/*
        abundarr[i][numA-1] = (float)(He[i]);
        abundarr[i][numA-1] = (float)((double)1.0 - He[i]);
*/
        i++;
    } while(!feof(afp));

    free(He);

    Abins = i;

    fclose(afp);
    afp = NULL;

    fp = fopen("abundances2", "w");

    for(i=0;i<Abins;i++) {
        for(j=0;j<numA;j++) {
            fprintf(fp,"%E ",(double)abundarr[i][j]);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
    fp = NULL;

    fp = fopen("c16r4_abun.dat", "w");
    if(!fp) printf("error opening file 'rung3g_50_abun.dat'\n");

    fprintf(fp,"%d\n",numA);
    fprintf(fp,"%d\n",Abins);
    for (j=0; j < numA; j++) fprintf(fp,"%d ",nparr[j]);
    fprintf(fp,"\n");
    for (j=0; j < numA; j++) fprintf(fp,"%d ",nnarr[j]);
    fprintf(fp,"\n");

    fprintf(fp,"%s ","r");
    for (j=0; j < numA; j++) fprintf(fp,"%s%-d ",elements[ nparr[j] ], nparr[j]+nnarr[j]);
    fprintf(fp,"\n");



    for(i=0;i<Abins;i++) {
        fprintf(fp,"%10E ",radbin[i]);
        for(j=0;j<numA;j++) {
            fprintf(fp,"%10E ", (double)abundarr[i][j]);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
    fp = NULL;

    /* read in the list of isotopes in the network */
    fp = fopen("network.isotopes", "r");
    if (!fp) printf("error opening files network.isotopes\n");

    fscanf(fp, "%d", &NNW);
    fscanf(fp, "%s\t%s", tmpchr,tmpchr);

    /* this way nwlist is indexed the same way as want */
    nwlist = (int **) malloc ( 2 * sizeof(int *) );
    nwlist[0] = (int *) malloc ( NNW * sizeof(int) );
    nwlist[1] = (int *) malloc ( NNW * sizeof(int) );

    for( i = 0; i < NNW; i++) {
        fscanf(fp, "%d", &nwlist[0][i]); /*protons*/
        fscanf(fp, "%d", &nwlist[1][i]); /*neutrons*/
    }

    fclose(fp);
    fp = NULL;

    fp = fopen("c16r4_abun_short.dat", "w");
    if(!fp) printf("error opening file 'run3g_50_abun_short.dat'\n");

    fprintf(fp,"%d\n",NNW);
    fprintf(fp,"%d\n",Abins);
    for (j=0; j < NNW; j++) fprintf(fp,"%d ",nwlist[0][j]);
    fprintf(fp,"\n");
    for (j=0; j < NNW; j++) fprintf(fp,"%d ",nwlist[1][j]);
    fprintf(fp,"\n");

    fprintf(fp,"%s ","r");
    for (j=0, m=0, n=0; j < numA; j++) {
        if(nparr[j] == nwlist[0][m] && nnarr[j] == nwlist[1][m]) {
            fprintf(fp,"%s%-d ",elements[ nparr[j] ], nparr[j]+nnarr[j]);
            m++;
        }
    }
    fprintf(fp,"\n");

    for(i=0;i<Abins;i++) {
        fprintf(fp,"%10E ",radbin[i]);
        for (j=0, m=0, n=0; j < numA; j++) {
            if(nparr[j] == nwlist[0][m] && nnarr[j] == nwlist[1][m]) {
                fprintf(fp,"%10E ", (double)abundarr[i][j]);
                m++;
            }
        }
        fprintf(fp,"\n");
    }

    return 0;
}
