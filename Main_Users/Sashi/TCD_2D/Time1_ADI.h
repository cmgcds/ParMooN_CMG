// ======================================================================
// instationary problem
// ======================================================================
#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Exzmple: Time1.h" << endl);
}

#define __ADI__

// exact solution
void Exact(double x, double y, double *values)
{
  double r2, t;

  t = TDatabase::TimeDB->CURRENTTIME;
  r2 = (x-0.5-0.2*t)*(x-0.5-0.2*t) + (y-0.5-0.1*t)*(y-0.5-0.1*t);

  if( r2 < 0.0625 )
    values[0] = 0.25 - sqrt(r2);
  else
    values[0] = 0;

  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_Space(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_Internal(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue_Space(int BdComp, double Param, double &value)
{
  value = 0;
}

// value of boundary condition
void BoundValue_Internal(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs_Internal(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  double a=0.2, b=0.1, c=0;

  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->RE_NR==0)
    eps =0.;
  else
    eps=1./TDatabase::ParamDB->RE_NR;



  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x= X[i];
//     y = Y[i];

    coeff[0] = eps; // internal diffusion
    coeff[1] = b; // internal convection
    coeff[2] = c; // internal reaction

    coeff[3] = 0;
  }
}

void BilinearCoeffs_Space(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  double a=0.2, b=0.1, c=0;

  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->RE_NR==0)
    eps =0.;
  else
    eps=1./TDatabase::ParamDB->RE_NR;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;  // space diffusion
    coeff[1] = a;    // space convection
    coeff[2] = c;    // space reaction

    coeff[3] = 0;     // space rhs (no need, since internal values are the rhs for space)
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double r2;

  r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);

  if( r2 < 0.0625 )
    values[0] = 0.25 - sqrt(r2);
  else
    values[0] = 0;
}


 void Generate1DUniformMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j;
  int *Lines;
  double len, h, x, y;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  len = End-Start;
  h = len/double(N_Cells);

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_Cells+1]; 

  x=Start;
  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(x, y);
    x =double(i+1)*h;
   }

  x=End;
  Vetrex[N_Cells] = new TVertex(x, y);

  CellTree = new TBaseCell*[N_Cells];


   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);

    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
//     ((TBaseCell *) SurfCellTree[i])->SetBd_Part(Bd_Part[i] );
   }

   Domain->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}


void FindinternalAdvection(int N, TBaseCell *Cell, TFEFunction1D **u, int N_InternalLevel, double *G,
                           double *QuadPtsRhs, double *QuadPtsRhsT)
{
  int  i, j, k, l, N_Points, N_BaseFunct, N_Sets=1;
  int *GlobalNumbers, *BeginIndex, *DOF;

  double *Weights, *zeta, X[20], AbsDetjk[20];
  double **origvaluesD0, **origvaluesD1, *sol, val;
  double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1;

  bool Needs2ndDer[1];

  TFESpace1D *FeSpace;
  FE1D FEId;
  TFE1D *Element;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Needs2ndDer[0] = FALSE;

  // assume that all internal levels use same FE Space
  FeSpace = u[0]->GetFESpace1D();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();

  FEId = FeSpace->GetFE1D(N, Cell);
  Element = TFEDatabase2D::GetFE1D(FEId);
  N_BaseFunct = Element->GetN_DOF();
  BaseFunct_ID = Element->GetBaseFunct1D_ID();

  bf = Element->GetBaseFunct1D();
  l = bf->GetPolynomialDegree();
  LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
  qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  qf1->GetFormulaData(N_Points, Weights, zeta);

  F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
  ((TLineAffin *)F_K)->SetCell(Cell);
  //((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);


  BaseFunct[0] = BaseFunct_ID;
  ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

  origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
  origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

  DOF = GlobalNumbers + BeginIndex[N];

  for(i=0; i<N_InternalLevel; i++)
   {
    sol = u[i]->GetValues();

    for(j=0; j<N_Points; j++)
     {
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      for(k=0;k<N_BaseFunct;k++)
       {
        val = orgD0[k] * sol[ DOF[k] ]; // for M
       }

      QuadPtsRhs[i*N_Points + j] = val;

     } //  for(j=0; 
   } // for(i=0*/

//   BilinearCoeffs();
  for(i=0; i<N_Points; i++)
   {
    G[i] = 0.2;
   } //  for(i=0; 


}





