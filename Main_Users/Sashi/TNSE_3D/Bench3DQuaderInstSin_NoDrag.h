// Benchmark problem, Quader, stationary

#define  __BENCH_SIN_NO_DRAG__
#define  __DOWNWIND__

void ExampleFile()
{
  OutPut("Example: Bench3DQuaderInstSin_NoDrag.h" << endl);
}

// ========================================================================
// initial condition
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
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
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  if( (fabs(x)<1e-6) || (fabs(x-2.5)<1e-6) )
  {
    value = 16*2.25*y*(0.41-y)*z*(0.41-z)*sin(Pi*t/8.0)/(0.41*0.41*0.41*0.41);
  }
  else
  {
    value = 0;
  }

}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
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
                           cell, SecondDer,
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

  cd *= -20/0.41;
  cl *= -20/0.41;

//  delete Derivatives[0];
  delete aux;
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
//      OutPut("cell " << i << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);
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
               // sort for x and y (first variant), cells in the center of the channel come first
               if (sx[j+1] - sx[j] < -1e-8)
                  change_now = 1;
               if ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(sy[j+1]-0.205) < fabs(sy[j]-0.205)-1e-8))
                  change_now = 1;
               if ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(fabs(sy[j+1]-0.205) - fabs(sy[j]-0.205))<1e-8)
                   &&  (fabs(sz[j+1]-0.205) < fabs(sz[j]-0.205)-1e-8))
                  change_now = 1;
               break;
            case 3:
               // sort for x and y (second variant)
               if ((sx[j+1] - sx[j] < -1e-8)||
                   ((fabs(sx[j+1]-sx[j]) < 1e-8) && (fabs(sy[j+1]-0.205) > fabs(sy[j]-0.205) )))
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
/*
   int ii;
   //  OutPut("changes " << i << " " << changes << endl);
   for (i=0;i<N_Cells;i++)
   {
//      OutPut("sort " << downwind[i] << endl);
      ii = downwind[i];
      //Cell = Coll->GetCell(ii);
//      OutPut("sort " << i << " " << downwind[i] << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);

      OutPut("sort " << i << " " << downwind[i] << " " << sx[i] << " " << sy[i] << " " << sz[i] << endl);
   }
*/  
   delete sx;
   delete sy;
   delete sz;

   OutPut("downwind numbering finished"<< endl);
   
   return;
}

