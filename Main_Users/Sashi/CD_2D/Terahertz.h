// ======================================================================
// smooth solution problem
// ======================================================================
// #include <ConvDiff2D.h>
#define __TERAHERTZ__

void ExampleFile()
{
  OutPut("Example: terahertz.h, c "  << endl) ;
 
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0]=0;
  values[1]=0; 
  values[2]=0;
  values[3]=0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
  
    switch(i)
     {
      case 0 : 
         cond = DIRICHLET;  
       break;     
       
      case  1:
         cond = NEUMANN;
      break;      
      case 2 : 
         cond = NEUMANN;  
       break;
            
      case  3:
         cond = NEUMANN;
      break;  
      
       default:
            OutPut("Unknown BoundCondition" << endl);
            exit(4711);;
     }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double alpha,  char_L;
  double x, y, rhsfact;

 
  alpha = TDatabase::ParamDB->P1;
  char_L = TDatabase::ParamDB->P4;
  rhsfact =  TDatabase::ParamDB->P7;
   
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    x = X[i];
    y = Y[i];
    
    
    if(x<=1)
     {
      coeff[4] = rhsfact*exp(alpha*y*char_L);
      
//       if(y<-40)
//       cout << y << " coeff[4] " <<  coeff[4] << endl;
      
     }
    else
     {coeff[4] = 0; }

    
  }
}


