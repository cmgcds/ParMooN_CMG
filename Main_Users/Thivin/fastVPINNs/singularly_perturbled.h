// ======================================================================
// Sine problem
// ======================================================================

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>

extern "C"
{
 #include <gridgen.h>    
  
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}



void ExampleFile()
{
  OutPut("Example: Square_inverse.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 
    cond = DIRICHLET;
  
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;

    if(BdComp==1)
        value = 0;
    else
        value = 0.0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;
  
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
    {
        coeff = coeffs[i];
        //double *param = parameters[i];

        coeff[0] = TDatabase::ParamDB->PE_NR; //eps
        // coeff[0] = (cos(x[i]) + cos(y[i])) ;
        coeff[1] = 1; // bx
        coeff[2] = 0; // by
        coeff[3] = 0; // c
        coeff[4] = 1 ; // f
    }
}



// COmmented out by thivin, due to the external library issue on ubuntu 20.04 
// Enable it only when you want to use tetgen for meshing

// // Generate the mesh using triangle.h
void  TriaReMeshGen(TDomain *&Domain)
{
 
} // TriaReMeshGen

