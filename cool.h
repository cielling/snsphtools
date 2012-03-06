#ifndef NISO
#define NISO 20 	/* number of isotopes tracked */
#endif

#ifndef COOLING
#define COOLING
float **tablep; //array to hold cooling curve table values
float **ionfracp; //array to hold ionfraction table values
#else
extern float **tablep;
extern float **ionfracp;
#endif

#ifndef BURNING
#define BURNING
int NNW;	/* number of isotopes in network */
int **inNW;
int nparr[NISO], nnarr[NISO];
#else
extern int NNW;
extern int **inNW;
extern int nparr[NISO];
extern int nnarr[NISO];
#endif

double find_ne(float abundarr[],int nparr[],int nnarr[],double rho,double temp,int Gridpts,int Nel);
double calc_lcool1(float abundarr[],int nparr[],int nnarr[],double rho,double temp,int Gridpts,int Nel,int extrapolate);
double calc_lcool2(double temp,double uint,double dens,int numu,int numD);

double analytic_cool(double temp);

void locate(float xx[], long Nel, float x, long *j);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

void polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy);

void init_CoolTable(int *Gridpts, int *Nel);
