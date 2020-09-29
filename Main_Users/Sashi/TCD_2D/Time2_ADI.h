// ======================================================================
// instationary problem
// ======================================================================
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: Time2_ADI.h" << endl);
}

#define __ADI__

// exact solution
void Exact(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 1+2*Param;
    break;
    case 1:
      value = 3+3*Param*t;
    break;
    case 2:
      value = 1+3*t+2*(1-Param);
    break;
    case 3:
      value = 1+3*t*(1-Param);
    break;
  } // endswitch
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=1, b=2, c=1;
  int i;
  double *coeff, *param;
  double x, y;
  double t=TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = a;
    coeff[2] = b;
    coeff[3] = c;

    coeff[4] = c*(1+2*x+3*t*y)
              +2*a
              +3*t*b
              +3*y;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 0;
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












