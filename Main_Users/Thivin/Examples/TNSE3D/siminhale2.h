// channel with circular cross section

#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>
#include <IsoBoundFace.h>
#include <MacroCell.h>
#include <BdSphere.h>
#include <cmath>
// #include <tetgen.h>

void ExampleFile()
{
#ifdef _MPI
	int rank;
	MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

	if (rank == TDatabase::ParamDB->Par_P0)
#endif
	{
		OutPut("Example: Siminhale2.h " << endl);
	}
}

void InitialU1(double x, double y, double z, double *values)
{
	values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
	values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
	values[0] = 0;
}

void InitialP(double x, double y, double z, double *values)
{
	values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double z, double *values)
{
	values[0] = 0;
	values[1] = 0;
	values[2] = 0;
	values[3] = 0;
	values[4] = 0;
}

void ExactU2(double x, double y, double z, double *values)
{
	values[0] = 0;
	values[1] = 0;
	values[2] = 0;
	values[3] = 0;
	values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{

	values[0] = 0;
	values[1] = 0;
	values[2] = 0;
	values[3] = 0;
	values[4] = 0;
}

void ExactP(double x, double y, double z, double *values)
{
	static double eps = 1 / TDatabase::ParamDB->RE_NR;

	values[0] = 0;
	values[1] = 0;
	values[2] = 0;
	values[3] = 0;
	values[4] = 0;
}

// 2 - walll
//  1 - out
// 0 - inlet

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
	TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;

	if (CompID == 0 || CompID == 2)
		cond = DIRICHLET;

	else
	{
		cond = NEUMANN;
	}
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
		double radius = 0.485;
	double U_Mean = 1.0;
	if (CompID == 0)
	{
		if (sqrt(pow((x), 2) + pow((z), 2)) < radius) // inside Circle
		{
			float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			value = (r- 0.5)*0.04;
			// cout << " Value : " << value <<endl;
			// value =4.405;
		}
		else
			value = 0.0;
	}
	
	else
		value = 0;
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
	// double radius = 4.92; // radius = 10mm/0.02
	// double x_c = 53.3948;
	// double z_c = 153.27;
	// double U_Mean = 1.0;

	double radius = 0.485;
	double U_Mean = 1.0;
	
	if (CompID == 0)
	{
		if (sqrt(pow((x), 2) + pow((z), 2)) < radius) // inside Circle
		{
			float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			value = U_Mean + (r- 0.5)*0.04;
			// cout << " Value : " << value <<endl;
			// value =4.405;
		}
		else
			value = 0.0;
	}
	
	else
	{
		value = 0.0;
	}


}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
		double radius = 0.485;
	double U_Mean = 1.0;
	if (CompID == 0)
	{
		if (sqrt(pow((x), 2) + pow((z), 2)) < radius) // inside Circle
		{
			float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			value = (r- 0.5)*0.04;
			// cout << " Value : " << value <<endl;
			// value =4.405;
		}
		else
			value = 0.0;
	}
		
	else
		value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
			   double **parameters, double **coeffs)
{
	static double eps = 1. / TDatabase::ParamDB->RE_NR;
	int i;
	double *coeff, x, y, z;
	for (i = 0; i < n_points; i++)
	{
		coeff = coeffs[i];

		coeff[0] = eps;
		coeff[1] = 0;						 // f1
		coeff[2] = 0;						 // -49050000; // ;  // f2:
		coeff[3] = 1.0 / 51.5414; // 3.22761; // f3
	}
}
