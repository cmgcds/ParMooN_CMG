// ======================================================================
// CD problem 1D
// ======================================================================


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


void BoundVales(int N_Inputs, double *Inn, double *Out) 
{
 Out[0] = 0; // L_min
 Out[1] = 0; // L_max
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
 

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // ld = Coords[5][i]; 

    // diffusion
    coeff[0] = eps_L;
    
    // convection
    coeff[1] = 1;
    
    // reaction term
    coeff[2] = 0;
    
    //rhs  
    coeff[3] = 1;
    // coeff[3] =  -k*(exp(-k*t))*sin(Pi*ld)*cos(Pi*lv)*cos(Pi*la) 
    //              + Pi*(exp(-k*t))*cos(Pi*ld)*cos(Pi*lv)*cos(Pi*la);
    // coeff[3] =  -k*ekt*sin(Pi*ld)  + Pi*ekt*cos(Pi*ld) + eps_L*Pi*Pi*sin(Pi*ld);
    // coeff[3] =  -k*ekt*sin(Pi*ld)  + eps_L*Pi*Pi*sin(Pi*ld);
        // cout << "ld: " << ld << " "<< coeff[3]  <<endl;
  }
}

