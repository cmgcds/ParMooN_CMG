// ======================================================================
// smooth solution problem
// ======================================================================
// #include <ConvDiff2D.h>
#define _OPTTOMOGRAPHY_

void ExampleFile()
{
  OutPut("Example: DiffOptTomography.h, c "  << endl) ;
 
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
 cond = FREESURF; 
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
  double *coeff, *param, mu_a, x_orig, y_orig;
  double  char_L;
  double x, y, r, rhsfact, sigma;

 

//   char_L = TDatabase::ParamDB->P4;
  rhsfact =  TDatabase::ParamDB->P7;
  mu_a = TDatabase::ParamDB->P1;
  sigma = TDatabase::ParamDB->P8; //   bubly spread radius
//   sigma *=sigma*2.;
  
  
  
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = mu_a;

    x = X[i];
    y = Y[i];
//     x_orig = x + 0.969;
//     y_orig = y -1.49e-7;
    x_orig = x + 41.7;
    y_orig = y -6.4e-6;
//     x_orig = x ;
//     y_orig = y  ;
    r= sqrt(x_orig*x_orig + y_orig*y_orig);
    
    if(r<sigma)
     {
      coeff[4] = rhsfact*exp(-(x_orig*x_orig + y_orig*y_orig)/sigma);
      
//       coeff[4] = rhsfact*(1. - r/sigma);
      
//         coeff[4] =  rhsfact*(1. + cos(Pi*r/sigma) )/2.;
      
      
      if(TDatabase::ParamDB->P4< coeff[4]) TDatabase::ParamDB->P4 = coeff[4];
     }
    else
     {coeff[4] = 0; }

    
  }
}


void RobinInt(TSquareMatrix2D *A, double *rhs, BoundCondFunct2D *BoundaryCondition)
{
 int i, j, k, l, m;
 int N_BaseFunct, *N_BaseFuncts;
 int N_Cells, N_Vertices, N_Edges;   
 int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF; 
 int *KCol, *RowPtr, *JointDOF, N_DOF;
 int comp, N_LinePoints, N_Active, index1, index2;
 int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;  
  
 double t0, t1, normn, val, Coeff, r;
 double *LineWeights, *zeta, *ValuesA;
 double  X_B[100], Y_B[100];
 double **uref;
 double uorig[MaxN_BaseFunctions2D];
  
 TBaseCell *cell;
 TFEDesc2D *FeDesc;
 BaseFunct2D *BaseFuncts;
 TCollection *Coll;
 TJoint *joint;
 TBoundEdge *boundedge;
 TBoundComp *BoundComp;
 BoundCond Cond0, Cond1;
 FE2D FEId;
 TFE2D *ele;
 RefTrans2D RefTrans;
 TRefTrans2D *F_K;
 TFESpace2D *fespace;
 BF2DRefElements RefElement;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;  
  
 BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
 N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

 fespace = A->GetFESpace();
 Coll = fespace->GetCollection();
 N_Cells = Coll->GetN_Cells();

 BeginIndex = fespace->GetBeginIndex();
 GlobalNumbers = fespace->GetGlobalNumbers();
 N_Active =  fespace->GetActiveBound();
  
 RowPtr = A->GetRowPtr();
 KCol = A->GetKCol();  
 ValuesA = A->GetEntries();  
 
 Coeff = 1./TDatabase::ParamDB->P3;

 for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge || joint->GetType() == BoundaryEdge)
      {
        boundedge = (TBoundEdge *)joint;
        BoundComp = boundedge->GetBoundComp();
        boundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
     } // endfor j
     
     
    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {
     FEId = fespace->GetFE2D(i, cell);

     for(j=0;j<N_IsoJoints;j++)
     {
      //cout << "Cell " << i << " has free surface." << endl;
      IJoint = JointNumbers[j];

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);      
      
      
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;
      } // endswitch     
      
      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);
 
      for(k=0;k<N_LinePoints;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         uorig[l] = uref[k][l]; 
 
        F_K->GetTangent(IJoint, zeta[k], t0, t1);  // old line 
        normn = sqrt(t0*t0+t1*t1);
        r = LineWeights[k]*normn*Coeff;

         for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];    
           index2 = RowPtr[TestDOF+1];

           for(m=0;m<N_BaseFunct;m++)
            {
             AnsatzDOF = DOF[m];
//              cout << AnsatzDOF << " -- " << TestDOF << endl;
             index1 = RowPtr[TestDOF];
  
             if(index1+1 == index2) continue;
  
             while(KCol[index1] != AnsatzDOF) index1++;              

             val = r*uorig[l]*uorig[m];
             ValuesA[index1] += val;  
//               cout << "A: " << TestDOF << " ";
//               cout << AnsatzDOF << " " << val << endl;   
            } //     for(m=0;m<N_BaseFunct;m++) 
          } //  for(l=0;l<N_BaseFunct;l++)
       } // for(k=0;k<N_LinePoints;k++
     } // for(j=0;j<N_IsoJoints;j++)
      
    } //  if(N_IsoJoints > 0)
     
  } // for(i=0;i<N_Cells;i++)
  
  
}//RobinInt
  



