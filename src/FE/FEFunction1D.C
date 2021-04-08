// =======================================================================
// @(#)FEFunction1D.C
//
// Class:       TFEFunction1D
// Purpose:     a function from a finite element space in 1D
//
// Author:      Sashikumaar Ganesan (17.05.2007)
//
// History:     start of implementation 17.05.2007
//
// =======================================================================

#include <Database.h>
#include <FEDatabase2D.h>
#include <FEFunction1D.h>
#include <string.h>
#include <AllRefTrans.h>
#include <MooNMD_Io.h>

#include <BaseFunct1D.h>
#include <NodalFunctional1D.h>
#include <stdlib.h>
/** constructor with vector initialization */
TFEFunction1D::TFEFunction1D(TFESpace1D *fespace1D, char *name,
                             char *description, double *values, int length)
{
  FESpace1D = fespace1D;

  Name = strdup(name);

  Description = strdup(description);

  Values = values;

  Length = length;
}

/** calculate the interpolation of an exact function */
void TFEFunction1D::Interpolate(DoubleFunctND *Exact)
{
  int i, j, k, l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF, N_Coord = 1;

  double *xi, *eta;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4], coord[2];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  RefTrans1D RefTrans, *RefTransArray;

#ifdef __2D__
  N_Coord = 2;
#endif

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble * N_DOFs);
  for (i = 0; i < N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);

    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X,
#ifdef __2D__
                                       Y,
#endif
                                       AbsDetjk);

    for (j = 0; j < N_Points; j++)
    {
      FctVal[0] = double(cell->GetRegionID());
      coord[0] = X[j];
#ifdef __2D__
      coord[1] = Y[j];
#endif
      Exact(N_Coord, coord, FctVal);
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(PointValues, FunctionalValues);

    DOF = GlobalNumbers + BeginIndex[i];

    for (j = 0; j < N_LocalDOFs; j++)
      Values[DOF[j]] = FunctionalValues[j];

  } // for(i=0;i<N_Cells
}

// Interpolate given 1D + 1d function, PBE
/** calculate the interpolation of an exact function */
void TFEFunction1D::Interpolate(int ConstCoord, double x, DoubleFunct2D *Exact)
{
  int i, j, k, l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  for (i = 0; i < N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Y, AbsDetjk);

    for (j = 0; j < N_Points; j++)
    {
      if (ConstCoord == 0)
      {
        Exact(x, Y[j], FctVal);
        //cout << " x " << x << " y " <<Y[j] << endl;
      }
      else
      {
        Exact(Y[j], x, FctVal);
      }

      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(PointValues, FunctionalValues);

    DOF = GlobalNumbers + BeginIndex[i];

    for (j = 0; j < N_LocalDOFs; j++)
      Values[DOF[j]] = FunctionalValues[j];

  } // for(i=0; i<N_Cells; i++)
}

// Interpolate given 2D + 1d function, PBE
/** calculate the interpolation of an exact function */
void TFEFunction1D::Interpolate(double x, double y, DoubleFunct3D *Exact)
{
  int i, j, k, l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  for (i = 0; i < N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);

    for (j = 0; j < N_Points; j++)
    {
      Exact(x, y, Z[j], FctVal);
      //cout << " x " << x << " y " <<Z[j] << endl;
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(PointValues, FunctionalValues);
    DOF = GlobalNumbers + BeginIndex[i];

    for (j = 0; j < N_LocalDOFs; j++)
      Values[DOF[j]] = FunctionalValues[j];

  } // for(i=0; i<N_Cells; i++)
}

// ** Compute the Mean Value of the solution and the gradient per cell **//
// if order = 0 , then it computes the mean Value of the Solution 
// if order = 1 , then it computes the mean value of gradient of the Solution 
void TFEFunction1D::InterpolateMeanValuePerCell(TFEFunction1D* MeanSolutionFEFunction,MultiIndex1D order)
{
  // variables declaration
  double *xi, *eta;
  bool Needs2ndDer[1];
  double Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double Intl_L[MaxN_PointsForNodal1D];
  int N_Points;             // Number of Quadraure Points in the  given Formula
  double* Weights;          // Quadrature Weights for the Choosen quadrature Formula
  double* zeta;             // Quadrature points in reference Domain.
  double** origvaluesD0;    // Values of shape Functions for the entire Quadrature Points ( for selected quad formula ) - 2D Array
  double** origvaluesD1;    // Values of Shape 
  double Mult;             //  Quadrature wt * det Jk
  double *orgD0, *orgD1 ;  //  Values of shape Functions at a Quadrature Point   
  int* DOF;                // Global to Local DOF Mapping per cell


  // Variables to be used in the cell Looping
	TCollection *Coll;
	TBaseCell *cell;
	FE1D FEId;
	TFE1D *Element;
	TBaseFunct1D *bf;
	QuadFormula1D LineQuadFormula;
	TQuadFormula1D *qf1;
	TRefTrans1D *F_K;
	BaseFunct1D BaseFunct_ID, BaseFunct[1];
  Needs2ndDer[0] = FALSE;

  Coll = FESpace1D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int* BeginIndex = FESpace1D->GetBeginIndex();
  int* GlobalNumbers = FESpace1D->GetGlobalNumbers();


  // get the FE Function Arrays 
  double* mean_Array = MeanSolutionFEFunction->GetValues();
  double* Solution   = Values;


  for ( int i = 0 ; i < N_Cells ; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
		int N_BaseFunct = Element->GetN_DOF();
		BaseFunct_ID = Element->GetBaseFunct1D_ID();

    //double Measure of Cell ( Length )
    double area = cell->GetMeasure();


    int l = bf->GetPolynomialDegree();
		LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2 * l);
		qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
		qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
		((TLineAffin *)F_K)->SetCell(cell);
		((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);   // Get the Original Co-orinates for the cell from xi values

    BaseFunct[0] = BaseFunct_ID;
		((TLineAffin *)F_K)->GetOrigValues(1, BaseFunct, N_Points, zeta, LineQuadFormula, Needs2ndDer);

    origvaluesD0 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
		origvaluesD1 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    // Get the Local to Global DOF Mapping 
    DOF = GlobalNumbers + BeginIndex[i];

    for ( int quadPt = 0 ; quadPt < N_Points ; quadPt++)
    {
      Mult = Weights[quadPt] * AbsDetjk[quadPt];
      orgD0 = origvaluesD0[quadPt];
			orgD1 = origvaluesD1[quadPt];


      for ( int j = 0 ; j < N_BaseFunct ; j++)   // looping over all the Nodes in a cell
      { 
          int GlobDOF =  DOF[j];
          if(order == 0)
            mean_Array[i] +=  orgD0[j]*Values[GlobDOF]*Mult;
          else if (order == 1)
            mean_Array[i] +=  orgD1[j]*Values[GlobDOF]*Mult;
          else
          {
            cout << " ERROR : Non valid Multi Index Variable Value " <<endl;
            cout << " FUNCTION : InterpolateMeanValuePerCell() --  File : FEFunction1D.C" <<endl;
            cout << " Error Description : As per Curreent Implementation, The mean Value can only be computed for Solution and its gradient " <<
                      "So the function currently accepts the value as D0 and D1 as parameters" <<endl;
          }
          
      }

    }
    if( fabs(area) > 1e-8 )  mean_Array[i] /= area ;
  }

  
}

// Interpolate given 2D/3D + 1d function's at nodal interpolation points, PBE
/** calculate the interpolation of an exact function */
void TFEFunction1D::InterpolateNodalPts(int N_Coord, double *Coords, DoubleFunctND *Exact, double *val)
{
  int i, j, k, l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, disp, N_GlobNodalPts;
  int *DOF, *IndexArray, *NodalPtIndex;

  double *xi, *eta;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4], x, y, z;

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  if (N_Coord > 0)
  {
    x = Coords[0];
    y = Coords[1];
#ifdef __3D__
    z = Coords[2];
#endif
  }
  N_Coord++; // this 1D Coord

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  NodalPtIndex = FESpace1D->GetIntlPtIndexOfPts();
  N_GlobNodalPts = FESpace1D->GetN_RootNodalPts();
  IndexArray = new int[N_GlobNodalPts];
  memset(IndexArray, 0, SizeOfInt * N_GlobNodalPts);
  memset(val, 0, SizeOfDouble * N_GlobNodalPts);

  disp = 0;
  for (i = 0; i < N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);

    for (j = 0; j < N_Points; j++)
    {
      //     //       Exact(x, y, Z[j], FctVal);
      Coords[N_Coord - 1] = Z[j];
      Exact(N_Coord, Coords, FctVal);
      k = NodalPtIndex[disp + j];
      val[k] += FctVal[0];
      IndexArray[k]++;
    }

    disp += N_Points;
  } // for(i=0; i<N_Cells; i++)

  for (i = 0; i < N_GlobNodalPts; i++)
  {
    if (IndexArray[i] == 0)
    {
      cout << "Error in TFEFunction1D::InterpolateNodalPts: " << IndexArray[i] << endl;
      exit(0);
    }
    val[i] /= (double)IndexArray[i];
  }
}

/** convert current grid to vector-values FE function */
void TFEFunction1D::GridToData()
{
  int i, j, k, l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  for (i = 0; i < N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);

    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Y, AbsDetjk);

    for (j = 0; j < N_Points; j++)
      PointValues[j] = Y[j];

    nf->GetAllFunctionals(PointValues, FunctionalValues);

    DOF = GlobalNumbers + BeginIndex[i];

    for (j = 0; j < N_LocalDOFs; j++)
      Values[DOF[j]] = FunctionalValues[j];

  } // for(i=0; i<N_Cells; i++)
}

TFEFunction1D::~TFEFunction1D()
{
  if (Name)
    delete Name;
  if (Description)
    delete Description;
}
