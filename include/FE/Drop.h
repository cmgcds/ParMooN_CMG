#include <math.h>
#include <Constants.h>
#include <MooNMD_Io.h>

// adap=0 - without adaptation; adap=1 - with adaptation
// a - parameter for grid movement, a=10 corresponds near the uniform grid;
int SchemeA(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);
int SchemeA_ax(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeA_ax1(double *r, double *z, double *h_s, double *h_n, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeA_ax11(double *r, double *z, double *h_s, double *h_n, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeA_ax2(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeB(double *r, double *z, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL);
        
int SchemeT4(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL,
        double *F, double *dF);
int SchemeT4_ax(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL,
        double *F, double *dF, int adap, double *a, int up);    
int SchemeT4_axNL(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL, 
        double *F, double *dF, int adap, double* a);   
void Solver_3diag(int N, double *c, double *d, double *e, double *b);

// a generator of a nonuniform grid
double S(double a, double t);

double dSdt(double a, double t);

// a new parameter for the generator of a nonuniform grid
double A(double curv, double h);
