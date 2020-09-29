// time dependent Navier-Stokes problem 
// test file for SSMUM

void ExampleFile()
{
  OutPut("Example: SSMUM_00.h" << endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
    values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
    values[0] = 0;
}
void InitialP(double x, double y, double *values)
{
    values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}
void ExactP(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    double omega, angle, nx, ny, mp_x, mp_y, x, y, r, velo;

    value = 0;
    // at the slit
    if (BdComp>=4)
    {
	// angular velocity 
	omega =  2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
	// current angle
	angle = TDatabase::ParamDB->SSMUM_ANGLE;
		// center of rotation
	mp_x = TDatabase::ParamDB->SSMUM_MP_X;
	mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
	// initially left boundary of the slit
	if (BdComp==4)
	{
	    //OutPut(Param << " ");
	    // compute original position
	    x = 0.45;
	    y = 0.35 + Param * 0.3;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity 
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (-1,0)
	    // is now rotated by angle
	    nx = -cos(angle);
	    ny = sin(angle);
	    if (Param>=0.5)
	    {
		value = -nx*velo;
	    }
	    else
	    {
		value = nx*velo;
	    }
	}
	// initially upper boundary
	if (BdComp==5)
	{
	    // compute original position
	    x = 0.45+ Param * 0.1;
	    y = 0.65;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (0,1)
	    // is now rotated by angle
	    nx = sin(angle);
	    ny = cos(angle);
	    value = ny * velo;
	}
	// initially right boundary of the slit
	if (BdComp==6)
	{
	    // compute original position
	    x = 0.55;
	    y = 0.65 - Param * 0.3;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (1,0)
	    // is now rotated by angle
	    nx = cos(angle);
	    ny = -sin(angle);
	    if (Param<=0.5)
	    {
		value = nx*velo;
	    }
	    else
	    {
		value = -nx*velo;
	    }		
	}
	// initially lower boundary
	if (BdComp==7)
	{
	    // compute original position
	    x = 0.55 - Param * 0.1;
	    y = 0.35;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (0,1)
	    // is now rotated by angle
	    nx = -sin(angle);
	    ny = -cos(angle);
	    value = ny * velo;
	}

    }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    double omega, angle, nx, ny, mp_x, mp_y, x, y, r, velo;
    value = 0;
    // at the slit
    if (BdComp>=4)
    {
	// angular velocity 
	omega =  2 * Pi * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
	// current angle
	angle = TDatabase::ParamDB->SSMUM_ANGLE;
	// center of rotation
	mp_x = TDatabase::ParamDB->SSMUM_MP_X;
	mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
	// initially left boundary of the slit
	if (BdComp==4)
	{
	    // compute original position
	    x = 0.45;
	    y = 0.35 + Param * 0.3;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (-1,0)
	    // is now rotated by angle
	    nx = -cos(angle);
	    ny = sin(angle);
	    if (Param<=0.5)
	    {
		value = ny*velo;
	    }
	    else
	    {
		value = -ny*velo;
	    }
	}
	// initially upper boundary
	if (BdComp==5)
	{
	    // compute original position
	    x = 0.45+ Param * 0.1;
	    y = 0.65;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (0,1)
	    // is now rotated by angle
	    nx = sin(angle);
	    ny = cos(angle);
	    value = -nx * velo;
	}
	// initially right boundary of the slit
	if (BdComp==6)
	{
	    // compute original position
	    x = 0.55;
	    y = 0.65 - Param * 0.3;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (1,0)
	    // is now rotated by angle
	    nx = cos(angle);
	    ny = -sin(angle);
	    if (Param<=0.5)
	    {
		value = ny*velo;
	    }
	    else
	    {
		value = -ny*velo;
	    }		
	}
	// initially lower boundary
	if (BdComp==7)
	{
	    // compute original position
	    x = 0.55 - Param * 0.1;
	    y = 0.35;
	    // distance from center of rotation
	    r = sqrt((x-mp_x)*(x-mp_x) + (y-mp_y)*(y-mp_y));
	    // velocity (THIS IS AN APPROXIMATION)
 	    velo = r * omega;
	    // compute normal vector
	    // original normal vector is (0,1)
	    // is now rotated by angle
	    nx = -sin(angle);
	    ny = -cos(angle);
	    value = -nx * velo;
	}

    }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  double nu=1/TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = 0;
    coeff[2] = 0;
  }
}


void CheckDOFs(TCollection *coll, TFESpace2D *velo_space,
	       TFESpace2D *press_space)
{
    int i, j, N_Cells, N_Edges;
    int *GlobalNumbers_velo, *BeginIndex_velo, *DOF, *DOF1;
    int *GlobalNumbers_press, *BeginIndex_press;
    double sx, sy, r;
    double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
    double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;
    double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
    double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
    TBaseCell *cell;

    GlobalNumbers_velo = velo_space->GetGlobalNumbers();
    BeginIndex_velo = velo_space->GetBeginIndex();
    GlobalNumbers_press = press_space->GetGlobalNumbers();
    BeginIndex_press = press_space->GetBeginIndex();

    N_Cells = coll->GetN_Cells();
    
    for(i=0;i<N_Cells;i++)       
    {               
	// next cell
	cell = coll->GetCell(i);
  
	N_Edges = cell->GetN_Edges();

	// barycenter
	sx = sy = 0;
	for (j=0;j<N_Edges; j++)
	{
	    sx += cell->GetVertex(j)->GetX();
	    sy += cell->GetVertex(j)->GetY();
	}
	sx /= N_Edges;
	sy /= N_Edges;
	// check if mesh cell can be changed
	r = (sx - mp_x)*(sx - mp_x)+(sy - mp_y)*(sy - mp_y);
	r = sqrt(r);
	if (r>outer_radius)
	    continue;
	if (r<inner_radius)
	    continue;
	OutPut("dof cell " << i << " velo ");
	DOF =  GlobalNumbers_velo + BeginIndex_velo[i];
	DOF1 =  GlobalNumbers_velo + BeginIndex_velo[i+1];
	for (j=0;j<7;j++)
	    OutPut(DOF[j] << " ");
	OutPut("press ");
	DOF =  GlobalNumbers_press + BeginIndex_press[i];
	DOF1 =  GlobalNumbers_press + BeginIndex_press[i+1];
	for (j=0;j<3;j++)
	    OutPut(DOF[j] << " ");
	OutPut(endl);

    }

}
