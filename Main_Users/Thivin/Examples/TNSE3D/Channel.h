// Navier-Stokes problem, 3D Channel flow
//
// u(x,y) = unknown
// p(x,y) = unknown
#include<cmath>

void ExampleFile()
{
#ifdef _MPI
    int rank;
    MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

    if (rank == TDatabase::ParamDB->Par_P0)
#endif
    {
        OutPut("Example: Channel3D.h " << endl);
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
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{ // For any example with mixed boundary condition specially
    // with NEUMANN bd cond INTERNAL_PROJECT_PRESSURE should be set as 0
    //  and one with full DIRICHLET as in DrivenCavity3D.h
    // INTERNAL_PROJECT_PRESSURE should be set as 1
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;

    if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        if (CompID == 1)
            cond = NEUMANN;
        else if (CompID == 0)
            cond = DIRICHLET;
        else if (CompID == 2)
            cond = DIRICHLET;
    }
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
        double radius = 0.5;
        double meanValue = 1;
        y = y+ 0.5;
        z = z+ 0.5;
        if (CompID == 0)
        {
            value = y * (1 - y) * z * (1 - z) ;
            value = value/(pow(0.5,4));
            value = value*meanValue;
        //       value=1;
        }
            
        else if (CompID == 1)
            value = 0;
        else if (CompID == 2)
            value = 0;
        else 
        {
            cout << " Error on Boundary Value " <<endl;
            exit(0);
        }
  
}


// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
    value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
    value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
    static double eps = 1 / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff;

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        coeff[0] = eps;
        coeff[1] = 0; // f1
        coeff[2] = -1; // f2
        coeff[3] = 0; // f3
    }
}
