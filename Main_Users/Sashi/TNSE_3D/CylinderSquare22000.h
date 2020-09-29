// Benchmark problem

#include <stdlib.h>
#include <Database.h>

#define  __BENCH__
//#define  __DOWNWIND__

// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: CylinderSquare22000.h " << TDatabase::ParamDB->P7
         << " % noise (only for U1 !!!, see [JK05])");
  if (!TDatabase::ParamDB->FRICTION_TYPE)
    {
     OutPut("no slip b.c." << endl);
    }
  else
     OutPut("slip+friction b.c." << endl);
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 22000; 
  /*if (TDatabase::ParamDB->CELL_MEASURE!=2)
  {
      TDatabase::ParamDB->CELL_MEASURE = 2;
      OutPut("CELL_MEASURE changed to " << 
	     TDatabase::ParamDB->CELL_MEASURE << endl);
	     }*/
  if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 100)
      
  {
      TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = 101;
      OutPut("TURBULENT_VISCOSITY_TYPE changed to " <<
	     TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE << endl);
  }
  OutPut("CYLINDER_22000_YPLUS_SIDES: " << TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES << endl); 
  OutPut("CYLINDER_22000_YPLUS_FRONT: " << TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT << endl); 
  OutPut("CYLINDER_22000_YPLUS_BACK:  " << TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK << endl); 
}
// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0.0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
     cond  =  SLIP_FRICTION_PENETRATION_RESISTANCE;
     TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;

    // outflow boundary condition
    if (fabs(x-2.5)<1e-6)
    {
	cond = NEUMANN;
	TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    }
 
    // inflow boundary condition
    if (fabs(x)<1e-6)
      cond  = DIRICHLET;
    
/*    // boundary conditions at the column
    if ((fabs(x-0.45)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
      {
	cond  = DIRICHLET;
	//cout << "left ";
      }
    if ((fabs(x-0.55)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
      {
	cond  = DIRICHLET;
	//cout << "right ";
      }
    if ((fabs(y-0.65)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
      {
	cond  = DIRICHLET;
	//cout << "lower ";
      }
    if ((fabs(y-0.75)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
      {
	cond  = DIRICHLET;
	//cout << "upper ";
      }
*/   
    // no slip b.c. on the bottom, top, column
    if (TDatabase::ParamDB->FRICTION_TYPE==0)
    {
	TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 3;
	/*// bottom and top 
	   if (fabs(z)<1e-6)
	   cond  = DIRICHLET;
	   if (fabs(z-0.4)<1e-6)
	   cond  = DIRICHLET;*/
	// boundary conditions at the column
	if ((fabs(x-0.45)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "left ";
	}
	if ((fabs(x-0.55)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "right ";
	}
	if ((fabs(y-0.65)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "lower ";
	}
	if ((fabs(y-0.75)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "upper ";
	}
        // this is for the cicular cylinder
	if (fabs((x-0.5)*(x-0.5) +(y-0.7)*(y-0.7)-2.5e-3)<=1e-6)
	    cond  = DIRICHLET;
    }
    else
    {
	// slip with friction boundary conditions on all walls
	TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 3;
    }
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
   double noise = TDatabase::ParamDB->P7/100.0;
   double eps = 1e-4;
    //if( (fabs(x)<1e-6) || (fabs(x-2.5)<1e-6) )
  if((fabs(x)<1e-6))
  {
      //if ((fabs(z)>eps)&&(fabs(z-0.4)>eps)&&(fabs(y)>eps)&&(fabs(y-1.4)>eps))
      	  value = 1 + noise * ((double)rand()/RAND_MAX-0.5);
	  //else
	  //value = 1;
  }
  else
  {
    value = 0;
  }

}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
   double noise = TDatabase::ParamDB->P7/1000.0;
  if((fabs(x)<1e-6))
  {
    value = noise* ((double)rand()/RAND_MAX-0.5);
    value = 0;
  }
  else
  {
    value = 0;
  }  
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
   double noise = TDatabase::ParamDB->P7/1000.0;
    if((fabs(x)<1e-6))
  {
    value = noise * ((double)rand()/RAND_MAX-0.5);
    value = 0;
  }
  else
  {
    value = 0;
  }    
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
      
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}
/** calculate characteristic values */
void GetCdCl(TFEFunction3D *u1fct, TFEFunction3D *u2fct,
             TFEFunction3D *u3fct,
             TFEFunction3D *pfct,
	     TFEFunction3D *u1oldfct, TFEFunction3D *u2oldfct,
	     double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Points,N_Faces,comp,ed_nr;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D];
  double Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  int N_LocalUsedElements;
  FE3D LocalUsedElements[2], CurrentElement;
  int *DOF;
  double **OrigFEValues, *Orig;
  boolean SecondDer[2] = { FALSE, FALSE };
  double *u1, *u2, *u3, *p, *u1old, *u2old;
  TFESpace3D *USpace, *PSpace;
  int *UGlobalNumbers, *UBeginIndex;
  int *PGlobalNumbers, *PBeginIndex;
  int *N_BaseFunct, N_Cells;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double value, value1, value2, value3, value4;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValues3[MaxN_BaseFunctions3D];
  double FEFunctValues4[MaxN_BaseFunctions3D];
  double FEFunctValues5[MaxN_BaseFunctions3D];
  double FEFunctValues6[MaxN_BaseFunctions3D];
  int N_DerivativesU = 4;
  double *Derivatives[MaxN_QuadPoints_3D];
  MultiIndex3D NeededDerivatives[4] = { D000, D100, D010, D001 };
  TFEFunction3D *vfct;
  double *v, nu = 1/TDatabase::ParamDB->RE_NR;
  double *Der, *aux;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
  TFE3D *eleCell;
  FE3D FEEle;
  TFEDesc3D *FEDesc;
  int N_DOF_Circ, *DOF_Circ;
  char VString[] = "v";
  double dt = TDatabase::TimeDB->TIMESTEPLENGTH;
  if (dt < 1e-8)
  {
    OutPut("time step to small " << endl);
    exit(4711);
  }

  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  u3 = u3fct->GetValues();
  u1old = u1oldfct->GetValues();
  u2old = u2oldfct->GetValues();
  p = pfct->GetValues();

  USpace = u1fct->GetFESpace3D();
  PSpace = pfct->GetFESpace3D();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  aux = new double [MaxN_QuadPoints_3D*19];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*19;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*SizeOfDouble);
  vfct = new TFEFunction3D(USpace, VString, VString, v, N_);

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    for(j=0;j<N_Faces;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryFace))
//          (joint->GetType() == IsoBoundface)) // boundary edge 
      {
        
        boundface = (TBoundFace *)joint;  
        BoundComp = boundface->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if ((comp>=4)&&(comp<=7)) 
          {
            FEEle = USpace->GetFE3D(i,cell);   // finite element of cell
            eleCell =  TFEDatabase3D::GetFE3D(FEEle); 
            FEDesc = eleCell->GetFEDesc3D();   // fe descriptor
            N_DOF_Circ = FEDesc->GetN_JointDOF(); // # local dofs on joints
            DOF_Circ = FEDesc->GetJointDOF(j); // local dofs on joint j
            DOF = UGlobalNumbers + UBeginIndex[i]; // pointer to global dofs
            for (k=0;k<N_DOF_Circ;k++)         // set fe on circle to 1 
              v[DOF[DOF_Circ[k]]] = 1;
          }
      }      
    }
  }
  
  cd = 0;
  cl = 0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 2;
    LocalUsedElements[0] = USpace->GetFE3D(i, cell);
    LocalUsedElements[1] = PSpace->GetFE3D(i, cell);

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta,
                           weights, X, Y, Z, AbsDetjk);

    // calculate all needed values of p 
    CurrentElement = LocalUsedElements[1];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = PGlobalNumbers + PBeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = p[DOF[l]];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
        value += FEFunctValues[l] * Orig[l];

      Derivatives[j][0] = value;
    }

    // calculate all needed values of u1, u2, u3 
    CurrentElement = LocalUsedElements[0];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];
    DOF = UGlobalNumbers + UBeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues1[l] = u1[DOF[l]];
      FEFunctValues2[l] = u2[DOF[l]];
      FEFunctValues3[l] = u3[DOF[l]];
      FEFunctValues4[l] = v[DOF[l]];
      FEFunctValues5[l] = u1old[DOF[l]];
      FEFunctValues6[l] = u2old[DOF[l]];
    }

    for(k=0;k<N_DerivativesU;k++)
    {
      OrigFEValues = TFEDatabase3D::
        GetOrigElementValues(BaseFunct,NeededDerivatives[k]);

      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value1 = 0;
        value2 = 0;
        value3 = 0;
        value4 = 0;
        for(l=0;l<N_;l++)
        {
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
          value3 += FEFunctValues3[l] * Orig[l];
          value4 += FEFunctValues4[l] * Orig[l];
        } // endfor l
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+5] = value2;
        Derivatives[j][k+9] = value3;
        Derivatives[j][k+13] = value4;
	if (k==0)
	{
	  value1 = 0;
	  value2 = 0;
	  for(l=0;l<N_;l++)
	    {
	      value1 += FEFunctValues5[l] * Orig[l];
	      value2 += FEFunctValues6[l] * Orig[l];
	    } // endfor l
	  Derivatives[j][17] = value1;
	  Derivatives[j][18] = value2;	  
	} // end if
      } // endfor j
    } // endfor k

// calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      // Der[0] = p
      // Der[1] = u1
      // Der[2] = u1_x
      // Der[3] = u1_y
      // Der[4] = u1_z 

      // nu * (u1_x*v_x+ u1_y*v_y + u1_z*v_z), v= (v,0,0)
      value1  = (Der[1]-Der[17])*Der[13]/dt 
	+ nu*(Der[2]*Der[14]+Der[3]*Der[15]+Der[4]*Der[16]);
      // (u1 * u1_x + u2* u1_y + u3* u1_z) * (1,0,0)
      value1 += (Der[1]*Der[2]+Der[5]*Der[3]+Der[9]*Der[4])*Der[13];
      // pressure times divergence of test function (1,0,0)
      value1 -= Der[0]*Der[14];

      value2  = (Der[5]-Der[18])*Der[13]/dt  + 
	nu*(Der[6]*Der[14]+Der[7]*Der[15]+Der[8]*Der[16]);
      value2 += (Der[1]*Der[6]+Der[5]*Der[7]+Der[9]*Der[8])*Der[13];
      value2 -= Der[0]*Der[15];

      cd += AbsDetjk[j]*weights[j] * value1;
      cl += AbsDetjk[j]*weights[j] * value2;
    }

  } // endfor i

  cd *= -50;
  cl *= -50;

  delete Derivatives[0];
  delete vfct;
  delete v;
}

/** calculate characteristic values */
void GetCdClOld(TFEFunction3D *u1fct, TFEFunction3D *u2fct,
             TFEFunction3D *u3fct,
             TFEFunction3D *pfct,
             double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Points,N_Faces,comp,ed_nr;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D];
  double Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  int N_LocalUsedElements;
  FE3D LocalUsedElements[2], CurrentElement;
  int *DOF;
  double **OrigFEValues, *Orig;
  boolean SecondDer[2] = { FALSE, FALSE };
  double *u1, *u2, *u3, *p;
  TFESpace3D *USpace, *PSpace;
  int *UGlobalNumbers, *UBeginIndex;
  int *PGlobalNumbers, *PBeginIndex;
  int *N_BaseFunct, N_Cells;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double value, value1, value2, value3, value4;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValues3[MaxN_BaseFunctions3D];
  double FEFunctValues4[MaxN_BaseFunctions3D];
  int N_DerivativesU = 4;
  double *Derivatives[MaxN_QuadPoints_3D];
  MultiIndex3D NeededDerivatives[4] = { D000, D100, D010, D001 };
  TFEFunction3D *vfct;
  double *v, nu = 1/TDatabase::ParamDB->RE_NR;
  double *Der, *aux;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
  TFE3D *eleCell;
  FE3D FEEle;
  TFEDesc3D *FEDesc;
  int N_DOF_Circ, *DOF_Circ;
  char VString[] = "v";

  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  u3 = u3fct->GetValues();
  p = pfct->GetValues();

  USpace = u1fct->GetFESpace3D();
  PSpace = pfct->GetFESpace3D();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  aux = new double [MaxN_QuadPoints_3D*20];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*20;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*SizeOfDouble);
  vfct = new TFEFunction3D(USpace, VString, VString, v, N_);

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    for(j=0;j<N_Faces;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryFace))
//          (joint->GetType() == IsoBoundface)) // boundary edge 
      {
        
        boundface = (TBoundFace *)joint;  
        BoundComp = boundface->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if ((comp>=4)&&(comp<=7)) 
          {
            FEEle = USpace->GetFE3D(i,cell);   // finite element of cell
            eleCell =  TFEDatabase3D::GetFE3D(FEEle); 
            FEDesc = eleCell->GetFEDesc3D();   // fe descriptor
            N_DOF_Circ = FEDesc->GetN_JointDOF(); // # local dofs on joints
            DOF_Circ = FEDesc->GetJointDOF(j); // local dofs on joint j
            DOF = UGlobalNumbers + UBeginIndex[i]; // pointer to global dofs
            for (k=0;k<N_DOF_Circ;k++)         // set fe on circle to 1 
              v[DOF[DOF_Circ[k]]] = 1;
          }
      }      
    }
  }
  
  cd = 0;
  cl = 0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 2;
    LocalUsedElements[0] = USpace->GetFE3D(i, cell);
    LocalUsedElements[1] = PSpace->GetFE3D(i, cell);

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta,
                           weights, X, Y, Z, AbsDetjk);

    // calculate all needed values of p 
    CurrentElement = LocalUsedElements[1];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = PGlobalNumbers + PBeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = p[DOF[l]];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
        value += FEFunctValues[l] * Orig[l];

      Derivatives[j][0] = value;
    }

    // calculate all needed values of u1, u2, u3 
    CurrentElement = LocalUsedElements[0];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];
    DOF = UGlobalNumbers + UBeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues1[l] = u1[DOF[l]];
      FEFunctValues2[l] = u2[DOF[l]];
      FEFunctValues3[l] = u3[DOF[l]];
      FEFunctValues4[l] = v[DOF[l]];
    }

    for(k=0;k<N_DerivativesU;k++)
    {
      OrigFEValues = TFEDatabase3D::
        GetOrigElementValues(BaseFunct,NeededDerivatives[k]);

      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value1 = 0;
        value2 = 0;
        value3 = 0;
        value4 = 0;
        for(l=0;l<N_;l++)
        {
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
          value3 += FEFunctValues3[l] * Orig[l];
          value4 += FEFunctValues4[l] * Orig[l];
        } // endfor l
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+5] = value2;
        Derivatives[j][k+9] = value3;
        Derivatives[j][k+13] = value4;
      } // endfor j
    } // endfor k

// calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      // Der[0] = p
      // Der[1] = u1
      // Der[2] = u1_x
      // Der[3] = u1_y
      // Der[4] = u1_z 

      // nu * (u1_x*v_x+ u1_y*v_y + u1_z*v_z), v= (v,0,0)
      value1  = nu*(Der[2]*Der[14]+Der[3]*Der[15]+Der[4]*Der[16]);
      // (u1 * u1_x + u2* u1_y + u3* u1_z) * (1,0,0)
      value1 += (Der[1]*Der[2]+Der[5]*Der[3]+Der[9]*Der[4])*Der[13];
      // pressure times divergence of test function (1,0,0)
      value1 -= Der[0]*Der[14];

      value2  = nu*(Der[6]*Der[14]+Der[7]*Der[15]+Der[8]*Der[16]);
      value2 += (Der[1]*Der[6]+Der[5]*Der[7]+Der[9]*Der[8])*Der[13];
      value2 -= Der[0]*Der[15];

      cd += AbsDetjk[j]*weights[j] * value1;
      cl += AbsDetjk[j]*weights[j] * value2;
    }

  } // endfor i

  cd *= -50;
  cl *= -50;

  delete Derivatives[0];
  delete vfct;
  delete v;
}

void DownwindNumberingCells(TCollection *Coll, int *downwind)
{
   int i,j,N_Cells, N_V, changes, change_now, tmp;
   double x, y, z, *sx, *sy, *sz;
   TBaseCell *cell;

   OutPut("downwind numbering started"<< endl);
      
   N_Cells = Coll->GetN_Cells();

   // take the numbering which is given by the GEO file
   if (TDatabase::ParamDB->SC_DOWNWIND_TYPE == 0)
   {
      for(i=0;i<N_Cells;i++)
         downwind[i] = i;
      return;
   }

   sx = new double[N_Cells];
   sy = new double[N_Cells];
   sz = new double[N_Cells];

   // compute center of mesh cells and store it
   for(i=0;i<N_Cells;i++)
   {
      cell = Coll->GetCell(i);
      // initialze the downwind array
      downwind[i] = i;
      N_V = cell->GetN_Vertices();
      x = 0;
      y = 0;
      z = 0;
      for (j=0;j<N_V;j++)
      {
         x += cell->GetVertex(j)->GetX();
         y += cell->GetVertex(j)->GetY();
         z += cell->GetVertex(j)->GetZ();
      }
      sx[i] = x/N_V;
      sy[i] = y/N_V;
      sz[i] = z/N_V;
//      OutPut(i << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);
   }

   changes = 1;
   i = 0;
   // sort the mesh cells
   while (changes)
   {
      i++;
      changes = 0;
      for (j=0;j<N_Cells-1;j++)
      {
         change_now = 0;
         switch(TDatabase::ParamDB->SC_DOWNWIND_TYPE)
         {
            case 1:
               // sort for x, cell with smaller x-value of its center comes first
               if (sx[j+1] - sx[j] < -1e-8)
                  change_now = 1;
               break;
            case 2:
               // sort for x and y (first variant)
               if ((sx[j+1] - sx[j] < -1e-8)||
                   ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(sy[j+1]-0.7) - fabs(sy[j]-0.7) < 0.0)))
                  change_now = 1;
               break;
            case 3:
               // sort for x and y (second variant)
               if ((sx[j+1] - sx[j] < -1e-8)||
                   ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(sy[j+1]-0.7) - fabs(sy[j]-0.7) > 0.0)))
                  change_now = 1;
               break;
            case 4:
               // sort for x, only at the beginning of the channel
                if ((sx[j+1] - sx[j] < -1e-8)&&(sx[j+1] <= 0.5 ))
                   change_now = 1;
               break;
                  
              
            default :
            {
               OutPut(" no routine for SC_DOWNWIND_TYPE = " << TDatabase::ParamDB->SC_DOWNWIND_TYPE 
                      << " defined" << endl);
               exit(4711);
            }
         }
         // if order of the mesh cells should be changed, do it
         // downwind[j] gets number of the mesh cells which comes 
         // first 
         if (change_now)
         {
            tmp = downwind[j];
            downwind[j] = downwind[j+1];
            downwind[j+1] = tmp;
            x = sx[j];
            sx[j] = sx[j+1];
            sx[j+1] = x;
            y = sy[j];
            sy[j] = sy[j+1];
            sy[j+1] = y;
            z = sz[j];
            sz[j] = sz[j+1];
            sz[j+1] = z;
            changes++;
         }
      } // endfor j
//      OutPut("changes " << i << " " << changes << endl);
   } // end while

   //  OutPut("changes " << i << " " << changes << endl);
//   for (i=0;i<N_Cells;i++)
//   {
//      OutPut("sort " << downwind[i] << endl);
//      OutPut("sort " << i << " " << downwind[i] << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);
//   }
   
   delete sx;
   delete sy;
   delete sz;

   OutPut("downwind numbering finished"<< endl);
   return;
}

void SetNoPenetrationValues(TSquareMatrix3D **SQMatricesA, TFEVectFunct3D *u, double h_min)
{
  int i, j, j0, j1, rows;  
  int *RowPtr, *KCol;
  double *u1, *u2, *u3, *Entries, comp;
  TSquareMatrix3D *sqmatrixA;

  
  u1 = u->GetComponent(0)->GetValues();
  u2 = u->GetComponent(1)->GetValues();
  u3 = u->GetComponent(2)->GetValues();

  comp = TDatabase::ParamDB->PENETRATION_CONSTANT * h_min 
      * TDatabase::TimeDB->TIMESTEPLENGTH/8.0;

  OutPut("comp " << comp << endl);

/*  sqmatrixA     = SQMatricesA[0];
  RowPtr        = sqmatrixA->GetRowPtr();
  KCol          = sqmatrixA->GetKCol();
  Entries       = sqmatrixA->GetEntries();
  rows          = sqmatrixA->GetN_Rows();

  for (i=0; i<rows; i++)
  {
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for (j=j0;j<j1;j++)
    {
	if (fabs(Entries[j])>comp)
	{
	    OutPut(" large1 " << KCol[j]);
	    u1[ KCol[j]] = 0;
	}
    }
  }
*/
  sqmatrixA     = SQMatricesA[1];
  RowPtr        = sqmatrixA->GetRowPtr();
  KCol          = sqmatrixA->GetKCol();
  Entries       = sqmatrixA->GetEntries();
  rows          = sqmatrixA->GetN_Rows();

  for (i=0; i<rows; i++)
  {
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for (j=j0;j<j1;j++)
    {
	if (fabs(Entries[j])>comp)
	{
	    //OutPut(" large2 " << KCol[j]);
	    u2[ KCol[j]] = 0;
	}
    }
  }

  sqmatrixA     = SQMatricesA[2];
  RowPtr        = sqmatrixA->GetRowPtr();
  KCol          = sqmatrixA->GetKCol();
  Entries       = sqmatrixA->GetEntries();
  rows          = sqmatrixA->GetN_Rows();

  for (i=0; i<rows; i++)
  {
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for (j=j0;j<j1;j++)
    {
	if (fabs(Entries[j])>comp)
	{
	    //OutPut(" large3 " << KCol[j]);
	    u3[ KCol[j]] = 0;
	}
    }
  }
}

/** calculate characteristic values */
void ComputeFrictionVelocities(TCollection *Coll,
			       TFEFunction3D *u1, TFEFunction3D *u2,
			       double *velo_friction, int &count)
{
    int i, j, k, N_V, N_Cells, N_Faces, comp;
    double h, values[8], no[4], val[4], x, y, z, eps=1e-5;
    TBaseCell *cell;
    TJoint *joint;
    TBoundFace *boundface;
    TBoundComp3D *BoundComp;
    TVertex *vertex;
    val[0] = val[1] = val[2] = val[3] = 0;
    no[0] = no[1] = no[2] = no[3] = 0;
    
// ########################################################################
// loop over all cells
// ########################################################################
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    for(j=0;j<N_Faces;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryFace))
//          (joint->GetType() == IsoBoundface)) // boundary edge 
      {
        
        boundface = (TBoundFace *)joint;  
        BoundComp = boundface->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if ((comp>=4)&&(comp<=7)) 
	{
	    N_V = cell->GetN_Vertices();
	    // find vertex on the boundary
	    for (k=0;k<N_V;k++)
	    {
		cell->GetVertex(k)->GetCoords(x,y,z);
		if ((fabs(x-0.45)<eps)||(fabs(x-0.55)<eps)||
		    (fabs(y-0.65)<eps)||(fabs(y-0.75)<eps))
		{
		    switch(comp)
		    {
			// front boundary
			case 4:
			    u1->FindGradientLocal(cell,i,x,y,z,values);
			    val[0] -= values[1];
			    no[0] += 1;
			    break;
			// upper boundary
			case 5:
			    u2->FindGradientLocal(cell,i,x,y,z,values);
			    val[1] += values[2];
			    no[1] += 1;
			    break;
			// back boundary
			case 6:
			    u1->FindGradientLocal(cell,i,x,y,z,values);
			    val[2] += values[1];
			    no[2] += 1;
			    break;
			// lower boundary
			case 7:
			    u2->FindGradientLocal(cell,i,x,y,z,values);
			    val[3] -= values[2];
			    no[3] += 1;
			    break;
		    }
		}
	    }

	}
      }
    }
  }
  if (count==0)
  {
      for (i=0;i<4;i++)
    	  velo_friction[i] = val[i]/no[i];
  }
  else
  {
      for (i=0;i<4;i++)
	  velo_friction[i] = count *(velo_friction[i]) /(count+1)  + (val[i]/no[i])/ (count+1);
  }

  OutPut(TDatabase::TimeDB->CURRENTTIME << " fric velo " << velo_friction[0]
	 << " " << velo_friction[1] << " " << velo_friction[2] << " " 
	 << velo_friction[3] << endl);

  OutPut(TDatabase::TimeDB->CURRENTTIME << " y+ " << sqrt(fabs(velo_friction[0])*TDatabase::ParamDB->RE_NR) << " " 
	 << sqrt(fabs(velo_friction[1])*TDatabase::ParamDB->RE_NR) << " "
	 << sqrt(fabs(velo_friction[2])*TDatabase::ParamDB->RE_NR) << " "
	 << sqrt(fabs(velo_friction[3])*TDatabase::ParamDB->RE_NR) << endl);

  count += 1;
}
/** calculate characteristic values */
void PreparePressureAtCylinder(TCollection *Coll, double* &press_cyl,
			       int &press_nodes)			      
{
    const int max_entries = 257;
    int i, j, k, N_V, N_Cells, N_Faces, comp, nodesx, nodesy, count = 0;
    int found, count1;
    double sx, sy, sx0, sy0, sx1, sy1, x, y, z, eps=1e-5;
    double press_cyl_x_tmp[max_entries], press_cyl_y_tmp[max_entries];
    double *press_cyl_x, *press_cyl_y;
    TBaseCell *cell;
    TJoint *joint;
    TBoundFace *boundface;
    TBoundComp3D *BoundComp;
    TVertex *vertex;
    
// ########################################################################
// loop over all cells
// ########################################################################
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      N_Faces=cell->GetN_Faces();
      for(j=0;j<N_Faces;j++)              // loop over all edges of cell
      {
	  joint=cell->GetJoint(j);
	  if ((joint->GetType() == BoundaryFace))
//          (joint->GetType() == IsoBoundface)) // boundary edge 
	  {
	      
	      boundface = (TBoundFace *)joint;  
	      BoundComp = boundface->GetBoundComp();  // get boundary component
	      comp=BoundComp->GetID();              // boundary id 
	      if ((comp>=4)&&(comp<=7)) 
	      {
		  N_V = cell->GetN_Vertices();
		  sx0 = sy0 = sx1 = sy1 = 0;
		  nodesx = nodesy = 0;
		  // find vertices on the boundary
		  // compute barycenter of boundary faces
		  for (k=0;k<N_V;k++)
		  {
		      cell->GetVertex(k)->GetCoords(x,y,z);
		      if ((fabs(x-0.45)<eps)||(fabs(x-0.55)<eps))
		      {
			  sx0 += x;
			  sy0 += y;
			  nodesx++;
		      }
		      if ((fabs(y-0.65)<eps)||(fabs(y-0.75)<eps))
		      {
			  sx1 += x;
			  sy1 += y;
			  nodesy++;
		      }
		  }
		  if (nodesx==4)
		  {
		      sx = sx0/4;
		      sy = sy0/4;
		  }
		  if (nodesy==4)
		  {
		      sx = sx1/4;
		      sy = sy1/4;
		  }
		  found = 0;
		  for (k=0;k<=count;k++)
		  {
		      if ((fabs(press_cyl_x_tmp[k]-sx)<eps) &&
			  (fabs(press_cyl_y_tmp[k]-sy)<eps))
		      {
			  found++;
			  break;
		      }
		  }
		  // new entry
		  if (!found)
		  {
		      press_cyl_x_tmp[count] = sx;
		      press_cyl_y_tmp[count] = sy;
		      count++;
		      if (count >= max_entries)
		      {
			  OutPut("PreparePressureAtCylinder: max_entries too small " 
				 << max_entries << endl);
			  exit(4711);
		      }
		  }
	      }
	  }
      }
  }
  // allocate array and initialize
  press_nodes = count;
  press_cyl = new double[3*count];
  press_cyl_x = press_cyl;
  press_cyl_y = press_cyl + count;
  memset(press_cyl_y+count, 0, count*SizeOfDouble);
  
  count1 = 0;
  // sort
  for (i=0;i<press_nodes;i++)
  {
      // front
      if (fabs(press_cyl_x_tmp[i]-0.45) < eps)
      {
	  press_cyl_x[count1] = press_cyl_x_tmp[i];
	  press_cyl_y[count1] = press_cyl_y_tmp[i];
	  // correct order
	  for (j=0;j<count1;j++)
	  {
	      if (press_cyl_y[count1] < press_cyl_y[j])
	      {
		  for (k=count1-1;k>=j;k--)
		      press_cyl_y[k+1] = press_cyl_y[k];
		  press_cyl_y[j] = press_cyl_y_tmp[i];
		  break;
	      }
	  }
	  count1++;
      }
  }
  count = count1;
  // sort
  for (i=0;i<press_nodes;i++)
  {
      // left side
      if (fabs(press_cyl_y_tmp[i]-0.75) < eps)
      {
	  press_cyl_x[count1] = press_cyl_x_tmp[i];
	  press_cyl_y[count1] = press_cyl_y_tmp[i];
	  // correct order
	  for (j=count;j<count1;j++)
	  {
	      if (press_cyl_x[count1] < press_cyl_x[j])
	      {
		  for (k=count1-1;k>=j;k--)
		      press_cyl_x[k+1] = press_cyl_x[k];
		  press_cyl_x[j] = press_cyl_x_tmp[i];
		  break;
	      }
	  }
	  count1++;
      }
  }
  count = count1;
  // sort
  for (i=0;i<press_nodes;i++)
  {
      // back side
      if (fabs(press_cyl_x_tmp[i]-0.55) < eps)
      {
	  press_cyl_x[count1] = press_cyl_x_tmp[i];
	  press_cyl_y[count1] = press_cyl_y_tmp[i];
	  // correct order
	  for (j=count;j<count1;j++)
	  {
	      if (press_cyl_y[count1] > press_cyl_y[j])
	      {
		  for (k=count1-1;k>=j;k--)
		      press_cyl_y[k+1] = press_cyl_y[k];
		  press_cyl_y[j] = press_cyl_y_tmp[i];
		  break;
	      }
	  }
	  count1++;
      }
  }
  // sort
  count = count1;
  for (i=0;i<press_nodes;i++)
  {
      // ride side
      if (fabs(press_cyl_y_tmp[i]-0.65) < eps)
      {
	  press_cyl_x[count1] = press_cyl_x_tmp[i];
	  press_cyl_y[count1] = press_cyl_y_tmp[i];
	  // correct order
	  for (j=count;j<count1;j++)
	  {
	      if (press_cyl_x[count1] > press_cyl_x[j])
	      {
		  for (k=count1-1;k>=j;k--)
		      press_cyl_x[k+1] = press_cyl_x[k];
		  press_cyl_x[j] = press_cyl_x_tmp[i];
		  break;
	      }
	  }
	  count1++;
      }
  }
  //for (i=0;i<press_nodes;i++)
  //  OutPut(i << " press " << press_cyl_x[i] << " " << press_cyl_y[i] << endl);
}

void PressureAtCylinder(TCollection *Coll, TFEFunction3D *p,
			double *press_cyl,
			int press_nodes, int &count)			      
{
    int i, j, k, N_V, N_Cells, N_Faces, comp, nodesx, nodesy;
    int *press_cyl_no;
    double sx, sy, sz, sx0, sy0, sz0, sx1, sy1, sz1, x, y, z, eps=1e-5, val[4];
    double *press_cyl_x, *press_cyl_y, *press_cyl_val, *press_cyl_tmp;
    TBaseCell *cell;
    TJoint *joint;
    TBoundFace *boundface;
    TBoundComp3D *BoundComp;
    TVertex *vertex;
    
// ########################################################################
// loop over all cells
// ########################################################################
    N_Cells = Coll->GetN_Cells();
 
    press_cyl_x = press_cyl;
    press_cyl_y = press_cyl + press_nodes;
    press_cyl_val = press_cyl_y + press_nodes;
    press_cyl_tmp = new double[press_nodes];
    memset(press_cyl_tmp, 0, press_nodes*SizeOfDouble);
    press_cyl_no = new int[press_nodes];
    memset(press_cyl_no, 0, press_nodes*SizeOfInt);
        
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      N_Faces=cell->GetN_Faces();
      for(j=0;j<N_Faces;j++)              // loop over all edges of cell
      {
	  joint=cell->GetJoint(j);
	  if ((joint->GetType() == BoundaryFace))
//          (joint->GetType() == IsoBoundface)) // boundary edge 
	  {
	      
	      boundface = (TBoundFace *)joint;  
	      BoundComp = boundface->GetBoundComp();  // get boundary component
	      comp=BoundComp->GetID();              // boundary id 
	      if ((comp>=4)&&(comp<=7)) 
	      {
		  N_V = cell->GetN_Vertices();
		  sx0 = sy0 = sx1 = sy1 = sz0 = sz1 = 0;
		  nodesx = nodesy = 0;
		  // find vertices on the boundary
		  // compute barycenter of boundary faces
		  for (k=0;k<N_V;k++)
		  {
		      cell->GetVertex(k)->GetCoords(x,y,z);
		      if ((fabs(x-0.45)<eps)||(fabs(x-0.55)<eps))
		      {
			  sx0 += x;
			  sy0 += y;
			  sz0 += z;
			  nodesx++;
		      }
		      if ((fabs(y-0.65)<eps)||(fabs(y-0.75)<eps))
		      {
			  sx1 += x;
			  sy1 += y;
			  sz1 += z;
			  nodesy++;
		      }
		  }
		  if (nodesx==4)
		  {
		      sx = sx0/4;
		      sy = sy0/4;
		      sz = sz0/4;
		  }
		  if (nodesy==4)
		  {
		      sx = sx1/4;
		      sy = sy1/4;
		      sz = sz1/4;
		  }
		  for (k=0;k<=press_nodes;k++)
		  {
		      if ((fabs(press_cyl_x[k]-sx)<eps) &&
			  (fabs(press_cyl_y[k]-sy)<eps))
		      {
			  p->FindGradientLocal(cell,i,sx,sy,sz,val);
			  press_cyl_tmp[k] += val[0];
			  press_cyl_no[k]++;
			  break;
		      }
		  }
	      }
	  }
      }
  }

  for (i=0;i<press_nodes;i++)
      press_cyl_tmp[i] /= press_cyl_no[i];

  for (i=0;i<press_nodes;i++)
  {
      press_cyl_val[i] = count *(press_cyl_val[i]) /(count+1)  
	  + press_cyl_tmp[i]/ (count+1);
      OutPut(TDatabase::TimeDB->CURRENTTIME << " p_cyl " << press_cyl_tmp[i] << " " 
	     << press_cyl_val[i] << endl);
  }
  count++;
  delete press_cyl_no;
  delete press_cyl_tmp;
}

void PrepareCenterlineVelocities(TCollection *Coll, double* &center_velo,
				 int &N_center_velo)
			      
{
    const int max_entries = 200;
    int i, j, k, N_Cells, N_V, comp, nodesx, nodesy, count = 0;
    int found, count1, found1;
    double x[8], y[8], z[8], eps=1e-5, xc;
    double center_velo_tmp[max_entries]; 
    TBaseCell *cell;
    
    for (i=0;i<max_entries;i++)
	center_velo_tmp[i] = -1;
    
// ########################################################################
// loop over all cells
// ########################################################################
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      found1 = 0;
      xc = 0;
      for (j=0;j<N_V;j++)
      {	  
	  cell->GetVertex(j)->GetCoords(x[j],y[j],z[j]);
	  if (fabs(y[j]-0.7)<eps)
	  {
	      found1 = 1;
	      found = 0;
	      xc += x[j];
	      for (k=0;k<=count;k++)
	      {
		  if ((fabs(center_velo_tmp[k]-x[j])<eps))
		  {
		      found++;
		      break;
		  }
	      }
	      // new entry
	      if (!found)
	      {
		  center_velo_tmp[count] = x[j];
		  count++;
		  if (count >= max_entries)
		  {
		      OutPut("PrepareCenterlineVelocities: max_entries too small " 
			     << max_entries << endl);
		      exit(4711);
		  }
	      }
	  }
      }
      // coordinate in the center
      if (found1)
      {
	  found = 0;
	  xc/=4;
	  for (k=0;k<=count;k++)
	  {
	      if ((fabs(center_velo_tmp[k]-xc)<eps))
	      {
		  found++;
		  break;
	      }
	  }
	  // new entry
	  if (!found)
	  {
	      center_velo_tmp[count] = xc;
	      count++;
	      if (count >= max_entries)
	      {
		  OutPut("PrepareCenterlineVelocities: max_entries too small " 
			 << max_entries << endl);
		  exit(4711);
	      }
	  }
      }
  }

  // allocate array and initialize
  N_center_velo = count;
  center_velo = new double[3*count];
  memset(center_velo+count, 0, 2*count*SizeOfDouble);
  
  count1 = 0;
  // sort
  for (i=0;i<N_center_velo;i++)
  {
      center_velo[count1] = center_velo_tmp[i];
      // correct order
      for (j=0;j<count1;j++)
      {
	  if (center_velo[count1] < center_velo[j])
	  {
	      for (k=count1-1;k>=j;k--)
		  center_velo[k+1] = center_velo[k];
	      center_velo[j] = center_velo_tmp[i];
	      break;
	  }
      }
      count1++;
  }  
}

void CenterlineVelocities(TCollection *Coll, TFEFunction3D *u1,
			  TFEFunction3D *u2, double *center_velo,
			  int N_center_velo, int &count_av)
			      
{
    int i, j, k, N_Cells, N_V, comp, nodesx, nodesy, count = 0;
    int found, count1, *center_velo_no, found1, j1;
    double x[8], y[8], z[8], eps=1e-5, val[5], xc;
    double *center_velo_x, *center_velo_y; 
    TBaseCell *cell;
    
    center_velo_x = new double[2*N_center_velo];
    memset(center_velo_x, 0, 2*N_center_velo*SizeOfDouble);
    center_velo_y = center_velo_x + N_center_velo;
    center_velo_no = new int[N_center_velo];
    memset(center_velo_no, 0, N_center_velo*SizeOfInt);
    
// ########################################################################
// loop over all cells
// ########################################################################
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      found1 = 0;
      for (j=0;j<N_V;j++)
      {	  
	  cell->GetVertex(j)->GetCoords(x[j],y[j],z[j]);
	  if (fabs(y[j]-0.7)<eps)
	  {
	      found1 = 1;
	      for (k=0;k<N_center_velo;k++)
	      {
		  if (fabs(center_velo[k]-x[j])<eps)
		  {
		      u1->FindGradientLocal(cell,i,x[j],y[j],z[j],val);
		      center_velo_x[k] += val[0];
		      u2->FindGradientLocal(cell,i,x[j],y[j],z[j],val+1);
		      center_velo_y[k] += val[1];	   
		      center_velo_no[k]++;
		      // top and bottom twice, only half of number of mesh cells
		      if ((fabs(z[j])<eps) || (fabs(z[j]- TDatabase::ParamDB->DRIFT_Z)<eps))
		      {
			  center_velo_x[k] += val[0];
			  center_velo_y[k] += val[1];	   
			  center_velo_no[k]++;
		      }
		  }
	      }
	  }
      }
      // coordinate in the center 
      if (found1)
      {
	  for (j=0;j<N_V;j++)
	  {
	      for (j1=j+1;j1<N_V;j1++)
	      {
		  if ((fabs(y[j]-0.7)<eps)&&(fabs(y[j1]-0.7)<eps)
		      && (fabs(z[j]-z[j1])<eps))
		  {
		      xc = (x[j]+x[j1])/2.0;
		      for (k=0;k<N_center_velo;k++)
		      {
			  if (fabs(center_velo[k]-xc)<eps)
			  {
			      u1->FindGradientLocal(cell,i,xc,y[j],z[j],val);
			      center_velo_x[k] += val[0];
			      u2->FindGradientLocal(cell,i,xc,y[j],z[j],val+1);
			      center_velo_y[k] += val[1];	   
			      center_velo_no[k]++;
			      // top and bottom twice, only half of number of mesh cells
			      if ((fabs(z[j])<eps) || (fabs(z[j]- TDatabase::ParamDB->DRIFT_Z)<eps))
			      {
				  center_velo_x[k] += val[0];
				  center_velo_y[k] += val[1];	   
				  center_velo_no[k]++;
			      }
			  }
		      }
		  }
	      }
	  }
      }
  }

  for (i=0;i<N_center_velo;i++)
  {
      center_velo_x[i]/=center_velo_no[i];
      center_velo_y[i]/=center_velo_no[i];
      center_velo[N_center_velo+i] = count_av * center_velo[N_center_velo+i]/(count_av+1)
	  + center_velo_x[i]/(count_av+1);
      center_velo[2*N_center_velo+i] = count_av * center_velo[2*N_center_velo+i]/(count_av+1)
	  + center_velo_y[i]/(count_av+1);
      OutPut(TDatabase::TimeDB->CURRENTTIME << " c_vel " << center_velo[i]  << " x " <<
	     center_velo_x[i] << " " << center_velo[N_center_velo+i] << " y " <<
	     center_velo_y[i] << " " << center_velo[2*N_center_velo+i] << endl);
  }
  count_av++;
  delete center_velo_x;
  delete center_velo_no;
}
void PrepareVelocityAtCylinder(TCollection *Coll, double* &center_velo,
			       int &N_center_velo)
			      
{
    const int max_entries = 200;
    int i, j, k, N_Cells, N_V, comp, nodesx, nodesy, count = 0;
    int found, count1, found1;
    double x[8], y[8], z[8], eps=1e-5, yc;
    double center_velo_tmp[max_entries]; 
    TBaseCell *cell;
    
    for (i=0;i<max_entries;i++)
	center_velo_tmp[i] = -1;
    
// ########################################################################
// loop over all cells
// ########################################################################
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      found1 = 0;
      yc = 0;
      for (j=0;j<N_V;j++)
      {	  
	  cell->GetVertex(j)->GetCoords(x[j],y[j],z[j]);
	  if (fabs(x[j]-0.5)<eps)
	  {
	      found1 = 1;
	      found = 0;
	      yc += y[j];
	      for (k=0;k<=count;k++)
	      {
		  if ((fabs(center_velo_tmp[k]-y[j])<eps))
		  {
		      found++;
		      break;
		  }
	      }
	      // new entry
	      if (!found)
	      {
		  center_velo_tmp[count] = y[j];
		  count++;
		  if (count >= max_entries)
		  {
		      OutPut("PrepareVelocityAtCylinder: max_entries too small " 
			     << max_entries << endl);
		      exit(4711);
		  }
	      }
	  }
      }
      // coordinate in the center
      if (found1)
      {
	  found = 0;
	  yc/=4;
	  for (k=0;k<=count;k++)
	  {
	      if ((fabs(center_velo_tmp[k]-yc)<eps))
	      {
		  found++;
		  break;
	      }
	  }
	  // new entry
	  if (!found)
	  {
	      center_velo_tmp[count] = yc;
	      count++;
	      if (count >= max_entries)
	      {
		  OutPut("PrepareVelocityAtCylinder: max_entries too small " 
			 << max_entries << endl);
		  exit(4711);
	      }
	  }
      }
  }

  // allocate array and initialize
  N_center_velo = count;
  center_velo = new double[2*count];
  memset(center_velo+count, 0, count*SizeOfDouble);
  
  count1 = 0;
  // sort
  for (i=0;i<N_center_velo;i++)
  {
      center_velo[count1] = center_velo_tmp[i];
      // correct order
      for (j=0;j<count1;j++)
      {
	  if (center_velo[count1] < center_velo[j])
	  {
	      for (k=count1-1;k>=j;k--)
		  center_velo[k+1] = center_velo[k];
	      center_velo[j] = center_velo_tmp[i];
	      break;
	  }
      }
      count1++;
  }  
}

void VelocityAtCylinder(TCollection *Coll, TFEFunction3D *u1,
			  double *center_velo,
			  int N_center_velo, int &count_av)
			      
{
    int i, j, k, N_Cells, N_V, comp, nodesx, nodesy, count = 0;
    int found, count1, *center_velo_no, found1, j1;
    double x[8], y[8], z[8], eps=1e-5, val[4], yc;
    double *center_velo_x;
    TBaseCell *cell;
    
    center_velo_x = new double[N_center_velo];
    memset(center_velo_x, 0, N_center_velo*SizeOfDouble);
    center_velo_no = new int[N_center_velo];
    memset(center_velo_no, 0, N_center_velo*SizeOfInt);
    
// ########################################################################
// loop over all cells
// ########################################################################
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
      cell = Coll->GetCell(i);
      N_V = cell->GetN_Vertices();
      found1 = 0;
      for (j=0;j<N_V;j++)
      {	  
	  cell->GetVertex(j)->GetCoords(x[j],y[j],z[j]);
	  if (fabs(x[j]-0.5)<eps)
	  {
	      found1 = 1;
	      for (k=0;k<N_center_velo;k++)
	      {
		  if (fabs(center_velo[k]-y[j])<eps)
		  {
		      u1->FindGradientLocal(cell,i,x[j],y[j],z[j],val);
		      center_velo_x[k] += val[0];
		      center_velo_no[k]++;
		      // top and bottom twice, only half of number of mesh cells
		      if ((fabs(z[j])<eps) || (fabs(z[j]- TDatabase::ParamDB->DRIFT_Z)<eps))
		      {
			  center_velo_x[k] += val[0];
			  center_velo_no[k]++;
		      }
		  }
	      }
	  }
      }
      // coordinate in the center 
      if (found1)
      {
	  for (j=0;j<N_V;j++)
	  {
	      for (j1=j+1;j1<N_V;j1++)
	      {
		  if ((fabs(x[j]-0.5)<eps)&&(fabs(x[j1]-0.5)<eps)
		      && (fabs(z[j]-z[j1])<eps))
		  {
		      yc = (y[j]+y[j1])/2.0;
		      for (k=0;k<N_center_velo;k++)
		      {
			  if (fabs(center_velo[k]-yc)<eps)
			  {
			      u1->FindGradientLocal(cell,i,x[j],yc,z[j],val);
			      center_velo_x[k] += val[0];
			      center_velo_no[k]++;
			      // top and bottom twice, only half of number of mesh cells
			      if ((fabs(z[j])<eps) || (fabs(z[j]- TDatabase::ParamDB->DRIFT_Z)<eps))
			      {
				  center_velo_x[k] += val[0];
				  center_velo_no[k]++;
			      }
			  }
		      }
		  }
	      }
	  }
      }
  }

  for (i=0;i<N_center_velo;i++)
  {
      center_velo_x[i]/=center_velo_no[i];
      center_velo[N_center_velo+i] = count_av * center_velo[N_center_velo+i]/(count_av+1)
	  + center_velo_x[i]/(count_av+1);
      OutPut(TDatabase::TimeDB->CURRENTTIME << " cyl_v " << center_velo[i]  << 
	     " " << center_velo_x[i] << " " << center_velo[N_center_velo+i] << endl);
  }
  count_av++;
  delete center_velo_x;
  delete center_velo_no;
}

