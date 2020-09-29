// flow over a wall mounted step

#define __CHANNELSTEP__

void ExampleFile()
{
  OutPut("Example: ChannelStepSlip3D.orig.h" << endl);
}

// ========================================================================
// exact solution
// ========================================================================
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
  // inflow 
  // if ((fabs(x)<1e-6)||(fabs(y)<1e-6) ||(fabs(y-10)<1e-6))
  if ((fabs(x)<1e-6))
    cond  = DIRICHLET;
  else
  {
    // outflow
    if (fabs(x-40)<1e-6)
    {
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    }
    else      
      // boundary conditions on the outer walls
    {
      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      //cond  = DIRICHLET;
    }
  }
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  //if( (fabs(x)<1e-6) || (fabs(x-2.5)<1e-6) )
  if((fabs(x)<1e-6))
  {
    value = 1;
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

/** calculate reattachment point */
void GetReattachmentLine(TFEFunction3D *u1fct,
                         double &reattachment_point)
{
  int i,j,k, N_Cells, l1, kk;
  double xi, eta, zeta;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  TFESpace3D *FESpace3D;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref, *u1;
  double *uzorig, *uzetaref;
  TJoint *joint;
  TBoundFace *boundedge;
  TBoundComp *BoundComp;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  
  int *Numbers, comp, found, N_Faces;
  double u[4], ux, uy, uz, x, y, z, *Values;
  double val, x_min, val_min,  x_max, val_max, eps=1e-12;
  double xv[4], valv[4];
  int *GlobalNumbers, *BeginIndex;

  FESpace3D = u1fct->GetFESpace3D();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  Values = u1fct->GetValues();  

  reattachment_point =0;

  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    found = 0;
    for(j=0;j<N_Faces;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if (joint->GetType() == BoundaryFace) // boundary face 
      {
        // compute point on the boundary face
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        
        for (l1=0;l1<TmpLen[j];l1++)
        {
          // get coordinates of the vertices belonging to the 
          // boundary face
          cell->GetVertex(TmpFV[j*MaxLen+l1])
            ->GetCoords(X[l1], Y[l1], Z[l1]);
        }
        // CORE DUMP MIT HPACC64, NUTZE HPACC20
        // face not in plane y=0
        if (!((fabs(Y[0])<1e-6)&&(fabs(Y[1])<1e-6)&&(fabs(Y[2])<1e-6)))
          continue;
        // face in front of the step
        if (X[0] < 6-1e-6)
          continue;
        
        // boundary face on plane y=0
        FE_ID = FESpace3D->GetFE3D(i, cell);
        FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
        RefTrans = FE_Obj->GetRefTransID();
        
        // get base function object
        bf = FE_Obj->GetBaseFunct3D();
        N_BaseFunct = bf->GetDimension();
        
        uorig = new double[N_BaseFunct];
        uxorig = new double[N_BaseFunct];
        uyorig = new double[N_BaseFunct];
        uzorig = new double[N_BaseFunct];
        
        uref = new double[N_BaseFunct];
        uxiref = new double[N_BaseFunct];
        uetaref = new double[N_BaseFunct];
        uzetaref = new double[N_BaseFunct];
 
        // set cell for reference transformation
        TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
   
        // compute tangential velocity for all vertices
        for (kk=0;kk<4;kk++)
        {
          TFEDatabase3D::GetRefFromOrig(RefTrans, X[kk], Y[kk], Z[kk], xi, eta, zeta);

          bf->GetDerivatives(D000, xi, eta, zeta, uref);
          bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
          bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
          bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);
      
          TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta, N_BaseFunct,
                                       uref, uxiref, uetaref, uzetaref, uorig, 
                                       uxorig, uyorig, uyorig);
      
          Numbers = GlobalNumbers + BeginIndex[i];
          
          u[kk] = 0;
          for(k=0;k<N_BaseFunct;k++)
          {
            val = Values[Numbers[k]];
            u[kk] += uorig[k]*val;
          }
          OutPut(" x " << X[kk] << " y " << Z[kk] << " u " << u[kk] << endl);
        } // endfor kk
        for (kk=0;kk<4;kk++)
        {
          for (k=kk+1;k<4;k++)
          {
            if (fabs(Z[kk]-Z[k])>1e-6)
              continue;
            // two vertices with the same y coordinates found
            // compute reattachment point on this line
            if ((X[kk] < X[k]) && (u[kk] <= 0 ) && (u[k] >=0) )
            {
              reattachment_point = X[kk] - u[kk]*(X[k]-X[kk])/(u[k]-u[kk]);
              OutPut("reattachment point : " << reattachment_point << " " << Z[kk] << endl);
            }
            if ((X[kk] > X[k]) && (u[kk] >= 0 ) && (u[k] <=0) )
            {
              reattachment_point = X[k] - u[k]*(X[kk]-X[k])/(u[kk]-u[k]);
              OutPut("reattachment point : " << reattachment_point << " " << Z[kk] << endl);
            }             
          } // endfor k
        } // endfor kk
      } // if BoundaryFace
    } 
    delete uorig;
    delete uxorig;
    delete uyorig;
    delete uzorig;
    delete uref;
    delete uxiref;
    delete uetaref;
    delete uzetaref;
  } // endfor i (cells)

}
