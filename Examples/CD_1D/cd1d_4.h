// ======================================================================
// CD problem 1D
// ======================================================================
#include<cmath>


void ExampleFile()
{
  OutPut("Example: CD1D.h" << endl) ;
}


// Exact Sol conditon
void Exact(int N_Coord, double *X, double *values)
{
 values[0] = 1;

  //  if(TDatabase::ParamDB->P2==1)
  // if(fabs(dist)>1)
  //  cout<< la << " : " << lv <<" : " << ld <<   " InitialPopulation  : " << InitialPopulation << endl;
}

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  // cond_Lmin = NEUMANN;
  // cond_Lmax = NEUMANN;

  cond_Lmin = DIRICHLET;
  cond_Lmax = DIRICHLET;
}

// BOUNDARY VALUES  
void BoundVales(int N_Inputs, double *Inn, double *Out) 
{
    double eps = 1.0/TDatabase::ParamDB->PE_NR;
 Out[0] = 1.0/(1.0 + exp(0.5/eps)); // L_min
 Out[1] = 1.0/(1.0 + exp(-0.5/eps)); // L_max
}

void BilinearCoeffs(int n_points, int N_Dim, double **Coords,
                    double **parameters, double **coeffs)
{
  int i;
  double eps_L, *coeff;
  // double inn[5], out[3];
 

  if(TDatabase::ParamDB->PE_NR)
    eps_L = 1.0/TDatabase::ParamDB->PE_NR;
  else
    eps_L = 0.;

    double *X = Coords[0];
  for(i=0;i<n_points;i++)
  {
    
    double x = X[i];
    // cout << " Val of X : " << x <<endl;
    double z = ( 0.5 - x) / eps_L;
    coeff = coeffs[i];
    // ld = Coords[5][i]; 

    // diffusion
    coeff[0] = eps_L;
    
    // convection
    coeff[1] = exp(-1.0/(x*x));
    
    // reaction term
    coeff[2] = 0;
    
    double ez = exp(z);
    double e2z = exp(2*z);

    double kk = (( (ez+1)*(ez) - 2*e2z )*exp(1.0/(x*x)) + (ez+1)* ez ) * exp(-1.0/(x*x));
    kk /= ( eps_L * pow((ez+1),3) );

    //rhs  
    coeff[3] = kk;
  }
}

