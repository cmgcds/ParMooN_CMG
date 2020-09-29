// ======================================================================
// instationary problem
// ======================================================================

#define __FUEL_CELL__
#define __STUTZEN__

/// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: brennStutzenCD_P.h" << endl); 
  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 1; // only rhs
}
// ======================================================================
// Sine problem 3D
// ======================================================================

// exact solution
void Exact(double x, double y, double z, double *values)
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
   if ((((x-17)*(x-17)+(y-120)*(y-120)-16 ) <=1e-6) && (fabs(z-0)<=1e-6))
     cond = DIRICHLET;
   else
     cond = NEUMANN;
}

// kind of boundary condition (for FE space needed)
void BoundConditionNSE(double x, double y, double z, BoundCond &cond)
{
   if ((((x-48)*(x-48)+(y-20)*(y-20)-16 ) <=1e-6) && (fabs(z-0)<=1e-6))
   {
      cond = NEUMANN;
   }
   else
      cond  = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  if ((((x-17)*(x-17)+(y-120)*(y-120)-16 ) <=1e-6) && (fabs(z-0)<=1e-6))
    value = 1;
  else
    value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
   if ((((x-17)*(x-17)+(y-120)*(y-120)-16 ) <=1e-6) && (fabs(z-0)<=1e-6))
     values[0] = 1;
   else
     values[0] = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = param[0];
    coeff[2] = param[1];
    coeff[3] = param[2];
    coeff[4] = 0;
    coeff[5] = 0;
  }
}

void GetOutFlowCells(TCollection *coll, int &N_OutFlowCells,
                     int* &OutFlowNumbers, int* &OutFlowFaces)
{
  int i,j,k;
  int N_Joints, N_Cells;
  TBaseCell *cell;
  TJoint *joint;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double sx, sy, sz, x, y, z;
  int IOutFlowCell;

  N_OutFlowCells = 0;
  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == BoundaryFace)
      {
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        sx = 0; sy = 0; sz = 0;
        for(k=0;k<TmpLen[j];k++)
        {
          cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(x, y, z);
          sx += x; sy += y; sz += z;
        } // endfor k
        sx /= TmpLen[j];
        sy /= TmpLen[j];
        sz /= TmpLen[j];
        if ((((sx-48)*(sx-48)+(sy-20)*(sy-20)-16 ) <=1e-6)
            && (fabs(sz-0)<=1e-6))
        {
          N_OutFlowCells++;
        } // endif
      } // endif
    } // endfor j
  } // endfor i

  OutPut("N_OutFlowCells: " << N_OutFlowCells << endl);

  OutFlowFaces = new int[N_OutFlowCells];
  OutFlowNumbers = new int[N_OutFlowCells];

  IOutFlowCell = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == BoundaryFace)
      {
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        sx = 0; sy = 0; sz = 0;
        for(k=0;k<TmpLen[j];k++)
        {
          cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(x, y, z);
          sx += x; sy += y; sz += z;
        } // endfor k
        sx /= TmpLen[j];
        sy /= TmpLen[j];
        sz /= TmpLen[j];
        if ((((sx-48)*(sx-48)+(sy-20)*(sy-20)-16 ) <=1e-6)
            && (fabs(sz-0)<=1e-6))
        {
          OutFlowNumbers[IOutFlowCell] = i;
          OutFlowFaces[IOutFlowCell] = j;
          IOutFlowCell++;
        } // endif
      } // endif
    } // endfor j
  } // endfor i
}

void GetOutFlow(int N_OutFlowCells, int *OutFlowCells,
                  int *OutFlowFaces, TFEFunction3D *conc,
                  double &volume, double &area)
{
  int i,j,k,l,m;
  int IJoint;
  TCollection* coll;
  TBaseCell *cell;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  QuadFormula2D FaceQuadFormula;
  TQuadFormula2D *qf2;
  int N_Points;
  double *FaceWeights, *t, *s;
  double **JointValues, *JointValue;
  int *N_BaseFuncts, N_BF;
  BaseFunct3D *BaseFuncts;
  double localvolume, localarea, v, a;
  int *BeginIndex, *GlobalNumbers, *DOF;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double *values;
  double X[4], Y[4], Z[4];
  double xc1, xc2, xc3, yc1, yc2, yc3, zc1, zc2, zc3;
  double nx, ny, nz;
  double detFK;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  values = conc->GetValues();

  fespace = conc->GetFESpace3D();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  coll = fespace->GetCollection();

  volume = 0;
  area = 0;

  for(i=0;i<N_OutFlowCells;i++)
  {
    cell = coll->GetCell(OutFlowCells[i]);
    IJoint = OutFlowFaces[i];
    
    CurrentElement = fespace->GetFE3D(OutFlowCells[i], cell);

    DOF = GlobalNumbers + BeginIndex[OutFlowCells[i]];

    cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    localvolume = 0;
    localarea = 0;

    for(l=0;l<TmpLen[IJoint];l++)
      cell->GetVertex(TmpFV[IJoint*MaxLen+l])->GetCoords(X[l], Y[l], Z[l]);
    
    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(CurrentElement);
    switch(TmpLen[IJoint])
    {
      case 3:
        // triangular face
        FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);
      break;

      case 4:
        // quadrilateral face
        FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
      break;
    }
    qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
    qf2->GetFormulaData(N_Points, FaceWeights, t, s);
    // generate data on reference mesh cell for the 2d face of 3d cell
    TFEDatabase3D::GetBaseFunct3DFromFE3D(CurrentElement)
      ->MakeRefElementData(FaceQuadFormula);
    // values of base functions in all quadrature points on face 
    JointValues = TFEDatabase3D::GetJointValues3D
      (BaseFuncts[CurrentElement], FaceQuadFormula, IJoint);

    N_BF = N_BaseFuncts[CurrentElement];

    switch(TmpLen[IJoint])
    {
      case 3:
        xc1 = X[1] - X[0];
        xc2 = X[2] - X[0];

        yc1 = Y[1] - Y[0];
        yc2 = Y[2] - Y[0];

        zc1 = Z[1] - Z[0];
        zc2 = Z[2] - Z[0];

        // normal vector
        nx = yc1*zc2 - zc1*yc2;
        ny = zc1*xc2 - xc1*zc2;
        nz = xc1*yc2 - yc1*xc2;
        // determinant of reference trafo
        detFK = sqrt(nx*nx + ny*ny + nz*nz);

        for(l=0;l<N_Points;l++)
        {
          JointValue = JointValues[l];

          v = 0;
          a = 0;
          for(m=0;m<N_BF;m++)
          {
            v += JointValue[m] * values[DOF[m]];
            a += JointValue[m];
          }

          localvolume += v*FaceWeights[l]*detFK;
          localarea += a*FaceWeights[l]*detFK;
        } // endfor l

      break;

      case 4:
        xc1=(-X[0] + X[1] + X[2] - X[3]) * 0.25;
        xc2=(-X[0] - X[1] + X[2] + X[3]) * 0.25;
        xc3=( X[0] - X[1] + X[2] - X[3]) * 0.25;

        yc1=(-Y[0] + Y[1] + Y[2] - Y[3]) * 0.25;
        yc2=(-Y[0] - Y[1] + Y[2] + Y[3]) * 0.25;
        yc3=( Y[0] - Y[1] + Y[2] - Y[3]) * 0.25;

        zc1=(-Z[0] + Z[1] + Z[2] - Z[3]) * 0.25;
        zc2=(-Z[0] - Z[1] + Z[2] + Z[3]) * 0.25;
        zc3=( Z[0] - Z[1] + Z[2] - Z[3]) * 0.25;

        for(l=0;l<N_Points;l++)
        {
          nx = (yc1+s[l]*yc3)*(zc2+t[l]*zc3)
              -(zc1+s[l]*zc3)*(yc2+t[l]*yc3);
          ny = (zc1+s[l]*zc3)*(xc2+t[l]*xc3)
              -(xc1+s[l]*xc3)*(zc2+t[l]*zc3);
          nz = (xc1+s[l]*xc3)*(yc2+t[l]*yc3)
              -(yc1+s[l]*yc3)*(xc2+t[l]*xc3);
          detFK = sqrt(nx*nx + ny*ny + nz*nz);

          JointValue = JointValues[l];

          v = 0;
          a = 0;
          for(m=0;m<N_BF;m++)
          {
            v += JointValue[m] * values[DOF[m]];
            a += JointValue[m];
          }

          localvolume += v*FaceWeights[l]*detFK;
          localarea += a*FaceWeights[l]*detFK;
        } // endfor l
      break;
    } // endswitch
  
    volume += localvolume;
    area += localarea;
  } // endfor i
}


