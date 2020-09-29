#include <constants.h>
#include <Enumerations.h>





void ALE_BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
	if(BdID == 3 || BdID == 0)
	{
		cond = DIRICHLET;
		//cout<<" BD ID for Dirichlet : "<<BdID<<endl;
	}
	else 
		cond = NEUMANN;
}

// Boundary Values for all the dimensions of the boundary domain. 
void ALE_BoundValue_X(int BdComp, double x, double y, double z, double &value)
{
    value = 0;
}

void ALE_BoundValue_Y(int BdComp, double x, double y, double z, double &value)
{
  		value = 0;
}

void ALE_BoundValue_Z(int BdComp, double x, double y, double z, double &value)
{
  value = 0;
}

void Grid_BoundCondition(int BdID, double x, double y, double z, BoundCond &cond)
{
	if(BdID == 3 || BdID == 0)
	{
		cond = DIRICHLET;
		//cout<<" BD ID for Dirichlet : "<<BdID<<endl;
	}
	else 
		cond = NEUMANN;
}


// Boundary Values for all the dimensions of the boundary domain. 
void Grid_BoundValue_X(int BdComp, double x, double y, double z, double &value)
{
    value = 0;
}

void Grid_BoundValue_Y(int BdComp, double x, double y, double z, double &value)
{
	if(BdComp == 3)
	{
		// value =  (5-x)*(5-x)  + (2.5-z)*(2.5-z) - 1.2*x*z ;
		// if(fabs(value) > 1e-6) value = value/48;
		// else value = 0;

		// value = 0.5 * sin( ( (3.141592654/2) + (x/10)  ) *3.1415926535) ;
		value 	=  x*0.15;	

	}
	else
	{
		value = 0;
	}
}

void Grid_BoundValue_Z(int BdComp, double x, double y, double z, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y,double *z,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;    coeff[1] = 0;  coeff[2] = 0; coeff[3] = 0;  coeff[4] = 0;
  }
}


void getFreeSurfaceBoundaryIds(std::vector<int>& Boundids)
{
	int size = 1;    // NOTE : Enter the Size of the Boundary ID in the given 
	// Boundids->resize(size);  // number of bd'ids to be considered as Free Surface
	Boundids.push_back(3);
}

void getFreeSlipBoundaryIds(std::vector<int>& FreeSlipBoundids)
{
	int size = 4;    // NOTE : Enter the Size of the Boundary ID in the given 
	FreeSlipBoundids.push_back(1);
	FreeSlipBoundids.push_back(2);
	FreeSlipBoundids.push_back(4);
	FreeSlipBoundids.push_back(5);
}

// THIVIN -
// This function "periodic boundary" is for storing the parameters that are needed for imposing the
// external boundary condition like tilting or swaying of the fluid

void externalBoundaryParameters(double& frequency, unsigned int& type, double& amplitude )
{
	 frequency = 1
	 ;   // ( cycles per second )
	
	// 1 is for Swaying ( planar Motion )
	// 2 is for tilting ( tilting motion )
	
	type = 1  ;

	 amplitude = 0.35;

}