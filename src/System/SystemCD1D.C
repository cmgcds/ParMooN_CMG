/** =======================================================================
* @class     TSystem1D
* @brief     stores the information of system CDR 1D
* @author    Sashikumaar Ganesan
* @date      12.12.2020
* @History 
* ======================================================================= */
#include <SystemCD1D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <DirectSolver.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <LineAffin.h>
#include <LinAlg.h>

#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <SquareMatrix.h>
#include <Matrix.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>

TSystemCD1D::TSystemCD1D(int N_L, double start, double end, BoundCond1D *boundConLminLMax, DoubleFunctND *BdValues, char *ParamFile) : TSystem1D(N_L, start, end, boundConLminLMax, BdValues, ParamFile)
{
	SqStructure = new TSquareStructure1D(FESpace1D);
	SqStructure->Sort();

	A_Intl = new TSquareMatrix1D(SqStructure);

	if (TDatabase::ParamDB->DISCTYPE == SUPG)
	{
		K_Intl = new TSquareMatrix1D(SqStructure);
	}

	// cout<<"test TSystemCD1D " <<endl;
	// exit(0);

} // SystemCD1D::TSystemCD1D

void TSystemCD1D::Init(CoeffFctND *bilinear)
{
	Bilinear = bilinear;
	char IString[] = "I";
	memset(Sol, 0, N_Dof * SizeOfDouble);
	FE_Function = new TFEFunction1D(FESpace1D, IString, IString, Sol, N_Dof);

	int N_DOF_Mean = this->GetN_Cells();
	// Allocate spaces for the Arrays
	meanSolution 		= new double[N_DOF_Mean]();
	meanSolutionGrad 	= new double[N_DOF_Mean]();

	cout << " here " <<endl;
	// Declare FESpace for Piecewise Const space P0
	FESpaceMean1D 			= new TFESpace1D(Coll_Intl, IString, IString, -1);
	
	// Declare FE Space Arrays 
	MeanSolution 			= new TFEFunction1D(FESpaceMean1D, IString, IString, meanSolution,		N_DOF_Mean);
	MeanSolutionGradient 	= new TFEFunction1D(FESpaceMean1D, IString, IString, meanSolutionGrad,	N_DOF_Mean);
} // Init



// To obtain the mean value of the solution and the gradient of the solution per each cell 

void TSystemCD1D::getMeanValueDerivatives()
{
	// Obtain the mean values
	FE_Function->InterpolateMeanValuePerCell(MeanSolution,D0);    			// Compute the mean value of the Solution per FE Cell
	FE_Function->InterpolateMeanValuePerCell(MeanSolutionGradient,D1);    // Compute the mean value of the Solution per FE Cell
	// cout << " mean Value Sol : " <<endl;
	// double* sol = MeanSolution->GetValues();
	// for ( int i = 0 ; i <  this->GetN_Cells(); i++)
	// {
	// 	cout << sol[i] << ",";
	// }
	// cout << endl;
	// cout<<endl;


	// cout << " mean Value Grad Sol : " <<endl;
	// double* sol1 = MeanSolutionGradient->GetValues();
	// for ( int i = 0 ; i <  this->GetN_Cells(); i++)
	// {
	// 	cout << sol1[i] << ",";
	// }
	// cout << endl;
	// cout<<endl;
	
}

//interpolation
void TSystemCD1D::Interpolate(DoubleFunctND *Exact)
{
	FE_Function->Interpolate(Exact);
} //TSystemCD1D::Interpolate

//solve the system
void TSystemCD1D::Solve()
{
	int N_Param = 0;
	double BDValue1, BDValue2;
	double *AdvectMatValues = A_Intl->GetEntries();
	if (cond_Lmin == DIRICHLET && TDatabase::ParamDB->DISCTYPE != DG)
		BDValue1 = Sol[0];

	if (cond_Lmax == DIRICHLET && TDatabase::ParamDB->DISCTYPE != DG)
		BDValue2 = Sol[N_Dof - 1];

	// init the matrices
	A_Intl->Reset();
	if (TDatabase::ParamDB->DISCTYPE == SUPG)
	{
		K_Intl->Reset();
	}

	// init rhs
	memset(Rhs, 0, N_Dof * SizeOfDouble);

	/** assemble matrices */
	switch (TDatabase::ParamDB->DISCTYPE)
	{
		//        case FD: // Finite difference, Il'in-Allen-Southwell scheme
		//            this->AssembleARhs_FD(G);
		//        break;
	case GALERKIN:
		this->AssembleARhs();
		break;
	case SDFEM: // SUPG
		this->AssembleARhs_SUPG();

		MatAdd(A_Intl, K_Intl, 1.0);
		break;

	case DG:
		//update B_nuc
		this->AssembleARhs_DG();
		break;
	default:
		Error("only FD, GALERKIN, SUPG, DG are implemented" << endl);
		Error("file: " << __FILE__ << " line " << __LINE__ << endl);
		exit(-1);
		break;
	}

	//Set BC
	this->SetDirichletBc();

	// //print matrix
	// int j, k, *RowPtr = A_Intl->GetRowPtr();
	// int *KCol = A_Intl->GetKCol();
	// double *ValuesA = A_Intl->GetEntries();
	//  for(j=0;j<N_Dof;j++)
	//   {
	//    for(k=RowPtr[j];k<RowPtr[j+1];k++)
	//     {
	//      cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
	//     }

	//      cout << "f: " << Rhs[j] <<endl;
	//   //  cout<<endl;
	//   }

	//solve the system
	DirectSolver(A_Intl, Rhs, Sol);

} //  TSystemCD1D::Solve

void TSystemCD1D::SetDirichletBc()
{
	int k;
	int *RowPtr, *KCol, begin, end;
	double *MatValues;

	RowPtr = A_Intl->GetRowPtr();
	KCol = A_Intl->GetKCol();
	MatValues = A_Intl->GetEntries();
	double BdVal[2];

	BundValues(0, NULL, BdVal);

	if (TDatabase::ParamDB->DISCTYPE != DG)
	{
		if (cond_Lmin == DIRICHLET)
		{
			begin = RowPtr[0];
			end = RowPtr[1];

			for (k = begin; k < end; k++)
			{
				if (KCol[k] == 0)
				{
					MatValues[k] = 1.;
				}
				else
				{
					MatValues[k] = 0.;
				}
			}

			//L_min value
			Rhs[0] = BdVal[0];
		}

		if (cond_Lmax == DIRICHLET)
		{
			begin = RowPtr[N_Dof - 1];
			end = RowPtr[N_Dof];

			for (k = begin; k < end; k++)
			{
				if (KCol[k] == N_Dof - 1)
				{
					MatValues[k] = 1.;
				}
				else
				{
					MatValues[k] = 0.;
				}
			}
			//L_max value
			Rhs[N_Dof - 1] = BdVal[1];
		}
	}
}

/**  Assembling A matrix Ilen-Southwell finite difference matrix*/
void TSystemCD1D::AssembleARhs_FD()
{
	int i, j, N;
	int begin, end, *RowPtr, *KCol;
	double *ValuesA, h, eps, q, a, b, c, beta;
	TCollection *Coll;
	TBaseCell *Cell;
	TVertex *Vetrex1, *Vetrex2;
	double g0 = 1;

	if (fabs(TDatabase::ParamDB->REACTOR_P3) > 0)
	{
		eps = 1.0 / TDatabase::ParamDB->REACTOR_P3;
	}
	else
	{
		eps = 1.e-15;
	}

	RowPtr = A_Intl->GetRowPtr();
	KCol = A_Intl->GetKCol();
	ValuesA = A_Intl->GetEntries();
	N = A_Intl->GetN_Columns();

	for (i = 0; i < N; i++)
	{
		if (i < (N - 1))
		{
			Cell = Coll->GetCell(i);
			Vetrex1 = Cell->GetVertex(0);
			Vetrex2 = Cell->GetVertex(1);
			h = fabs(Vetrex1->GetX() - Vetrex2->GetX());
		}

		// 2nd order Ilin-Allen Southwell scheme
		if (fabs(g0) > 0.)
		{
			q = g0 * h / (2. * eps);
			beta = eps * (q / tanh(q));
		}

		begin = RowPtr[i];
		end = RowPtr[i + 1];

		//BD
		a = -beta / (h * h) - g0 / (2 * h);
		b = 2 * beta / (h * h);
		c = -beta / (h * h) + g0 / (2 * h);

		for (j = begin; j < end; j++)
			if (i == (KCol[j] + 1)) //lower
			{
				ValuesA[j] = a;
			}
			else if (i == (KCol[j])) //diagonal
			{
				ValuesA[j] = b;
			}
			else if (i == (KCol[j] - 1)) //upper
			{
				ValuesA[j] = c;
			}
		//  cout<<endl;
	}

	// //print matrix
	//  for(j=0;j<N_Dof;j++)
	//   {
	//    begin = RowPtr[j];
	//    end = RowPtr[j+1];
	//    for(int k=begin;k<end;k++)
	//     {
	//      cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
	//     }

	//     //  cout << "f: " << rhs[j ] <<endl;
	//    cout<<endl;
	//   }
	//  exit(0);
}

// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TSystemCD1D::AssembleARhs()
{
	int i, j, k, l, N_Cells_Internal, N_BaseFunct;
	int N_Points, N_Sets = 1, *GlobalNumbers, *BeginIndex, *DOF;
	int TestDOF, begin, end, *RowPtr, *KCol;

	double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
	double LocMatrixA[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocRhs[MaxN_BaseFunctions1D];
	double **origvaluesD0, **origvaluesD1, Mult;
	double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
	double c0, c1, g0, rhsval, val, len = 0.;
	double **aux, **coeff, *Coeff;
	double **Coords;

	bool Needs2ndDer[1];

	TCollection *Coll;
	TBaseCell *Cell;
	FE1D FEId;
	TFE1D *Element;
	TBaseFunct1D *bf;
	QuadFormula1D LineQuadFormula;
	TQuadFormula1D *qf1;
	TRefTrans1D *F_K;
	BaseFunct1D BaseFunct_ID, BaseFunct[1];

	Coll = FESpace1D->GetCollection();
	GlobalNumbers = FESpace1D->GetGlobalNumbers();
	BeginIndex = FESpace1D->GetBeginIndex();

	RowPtr = A_Intl->GetRowPtr();
	KCol = A_Intl->GetKCol();
	ValuesA = A_Intl->GetEntries();

	N_Cells_Internal = Coll->GetN_Cells();
	Needs2ndDer[0] = FALSE;

	aux = new double *[MaxN_QuadPoints_1D];
	coeff = new double *[MaxN_QuadPoints_1D];
	Coords = new double *[N_Coord + 1];

	for (i = 0; i < N_Coord + 1; i++)
		Coords[i] = new double[MaxN_QuadPoints_1D];

	for (i = 0; i < MaxN_QuadPoints_1D; i++)
	{
		for (int ii = 0; ii < N_Coord; ++ii)
			Coords[ii][i] = 0.;

		aux[i] = new double[6];
		coeff[i] = new double[6];
	}

	Intl_L = Coords[N_Coord];

	for (i = 0; i < N_Cells_Internal; i++)
	{
		Cell = Coll->GetCell(i);
		FEId = FESpace1D->GetFE1D(i, Cell);
		Element = TFEDatabase2D::GetFE1D(FEId);
		bf = Element->GetBaseFunct1D();
		N_BaseFunct = Element->GetN_DOF();
		BaseFunct_ID = Element->GetBaseFunct1D_ID();

		l = bf->GetPolynomialDegree();
		LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2 * l);
		qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
		qf1->GetFormulaData(N_Points, Weights, zeta);

		F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
		((TLineAffin *)F_K)->SetCell(Cell);
		((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

		Bilinear(N_Points, N_Coord + 1, Coords, aux, coeff);

		BaseFunct[0] = BaseFunct_ID;
		((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta, LineQuadFormula, Needs2ndDer);

		origvaluesD0 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
		origvaluesD1 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

		memset(LocMatrixA, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocRhs, 0, N_BaseFunct * SizeOfDouble);

		DOF = GlobalNumbers + BeginIndex[i];

		for (j = 0; j < N_Points; j++)
		{
			Coeff = coeff[j];
			c0 = Coeff[0];	   // diffusion
			g0 = Coeff[1];	   // convection
			c1 = Coeff[2];	   // reaction term
			rhsval = Coeff[3]; //rhs

			// cout<< " c  " << c0 << " g0 " <<g0 << " c1 " << c1 <<endl;
			Mult = Weights[j] * AbsDetjk[j];
			orgD0 = origvaluesD0[j];
			orgD1 = origvaluesD1[j];

			len += Mult;
			for (k = 0; k < N_BaseFunct; k++)
			{
				test0 = orgD0[k];
				test1 = orgD1[k];
				// cout<< " test0  " << test0 << " orgD1  " <<test1 << " AbsDetjk " << AbsDetjk[j] <<endl;
				LocRhs[k] += Mult * rhsval * test0;

				for (l = 0; l < N_BaseFunct; l++)
				{
					ansatz0 = orgD0[l];
					ansatz1 = orgD1[l];

					// val  = g0*ansatz1*test0;// convective form
					val = -g0 * ansatz0 * test1; // conservative form
					val += c0 * ansatz1 * test1; // difusive term
					val += c1 * ansatz0 * test0;
					//  cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
					LocMatrixA[k * N_BaseFunct + l] += (Mult * val);
				}
			}
		}

		// add to global matrices
		for (j = 0; j < N_BaseFunct; j++)
		{
			TestDOF = DOF[j];
			Rhs[TestDOF] += LocRhs[j];

			begin = RowPtr[TestDOF];
			end = RowPtr[TestDOF + 1];
			for (k = begin; k < end; k++)
			{
				for (l = 0; l < N_BaseFunct; l++)
				{
					if (KCol[k] == DOF[l])
					{
						ValuesA[k] += LocMatrixA[j * N_BaseFunct + l];
						break;
					}
				} // for(m=0;m<N_BaseFunct_low
			}	  // for(n=begin;n<end;n++)

			//  //update boundary flux
			//  if(TestDOF==0 && cond_Lmin==NEUMANN)
			//   {
			//   for(k=begin;k<end;k++)
			//     {
			//      if(KCol[k] == 0 )
			//       { ValuesA[k] += -1; }
			//   //    else
			//   //     { ValuesM[k] = 0.; }
			//    }
			//   } //  if(cond_Lmin==DIRIC

		} // for(l=0;l<N_BaseFunct_low
	}	  // for(i=0; i<N_Cells_Internal; i++)

	//  cout<< " len " << len << endl;
	//  exit(0);
	//   //print matrix
	//    for(j=0;j<N_Dof;j++)
	//     {
	//      begin = RowPtr[j];
	//      end = RowPtr[j+1];
	//      for(k=begin;k<end;k++)
	//       {
	//        cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
	//       }

	//        cout << "f: " << Rhs[j ] <<endl;
	//     //  cout<<endl;
	//     }
	//  exit(0);
	for (i = 0; i < MaxN_QuadPoints_1D; i++)
	{
		delete[] aux[i];
		delete[] coeff[i];
	}

	delete[] aux;
	delete[] coeff;

	for (i = 0; i < N_Coord + 1; i++)
		delete[] Coords[i];

	delete[] Coords;

	//  printf("    AssembleARhs  \n" );
	//   MPI_Finalize();
	// exit(0);
} // void TSystemCD1D::AssembleARhs(d

/**  Assembling A matrix */
void TSystemCD1D::AssembleARhs_SUPG()
{
	int i, j, k, l, N_Cells_Internal, N_BaseFunct;
	int N_Points, N_Sets = 1, *GlobalNumbers, *BeginIndex, *DOF;
	int TestDOF, begin, end, *RowPtr, *KCol;

	double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
	double LocMatrixA[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocMatrixK[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocRhs[MaxN_BaseFunctions1D];
	double **origvaluesD0, **origvaluesD1, **origvaluesD2, Mult;
	double *orgD0, *orgD1, *orgD2, test0, test1, ansatz0, ansatz1, ansatz2, *ValuesA, *ValuesS, *ValuesK;
	double c0, c1, g0, rhsval, val, len = 0., bgradv, hE, beta, Pe_K;
	double **aux, **coeff, *Coeff, **Coords;

	bool Needs2ndDer[1];

	TCollection *Coll;
	TBaseCell *Cell;
	FE1D FEId;
	TFE1D *Element;
	TBaseFunct1D *bf;
	QuadFormula1D LineQuadFormula;
	TQuadFormula1D *qf1;
	TRefTrans1D *F_K;
	BaseFunct1D BaseFunct_ID, BaseFunct[1];
	TVertex *Vetrex1, *Vetrex2;

	double D_L = TDatabase::ParamDB->REACTOR_P3;
	double delta0 = TDatabase::ParamDB->DELTA0;
	double delta1 = TDatabase::ParamDB->DELTA1;
	int order = TDatabase::ParamDB->ANSATZ_ORDER;
	if (D_L < 1.e-12)
		D_L = 1.e-12;

	Coll = FESpace1D->GetCollection();
	GlobalNumbers = FESpace1D->GetGlobalNumbers();
	BeginIndex = FESpace1D->GetBeginIndex();

	RowPtr = A_Intl->GetRowPtr();
	KCol = A_Intl->GetKCol();
	ValuesA = A_Intl->GetEntries();
	ValuesK = K_Intl->GetEntries();

	N_Cells_Internal = Coll->GetN_Cells();
	Needs2ndDer[0] = TRUE;

	aux = new double *[MaxN_QuadPoints_1D];
	coeff = new double *[MaxN_QuadPoints_1D];

	Coords = new double *[N_Coord + 1];

	for (i = 0; i < N_Coord + 1; i++)
		Coords[i] = new double[MaxN_QuadPoints_1D];

	for (i = 0; i < MaxN_QuadPoints_1D; i++)
	{
		for (int ii = 0; ii < N_Coord; ++ii)
			Coords[ii][i] = 0;

		aux[i] = new double[6];
		coeff[i] = new double[6];
	}

	Intl_L = Coords[N_Coord];

	for (i = 0; i < N_Cells_Internal; i++)
	{
		Cell = Coll->GetCell(i);
		FEId = FESpace1D->GetFE1D(i, Cell);
		Element = TFEDatabase2D::GetFE1D(FEId);
		bf = Element->GetBaseFunct1D();
		N_BaseFunct = Element->GetN_DOF();
		BaseFunct_ID = Element->GetBaseFunct1D_ID();

		l = bf->GetPolynomialDegree();
		LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(20 * l);
		qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
		qf1->GetFormulaData(N_Points, Weights, zeta);

		F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
		((TLineAffin *)F_K)->SetCell(Cell);
		((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

		Bilinear(N_Points, N_Coord + 1, Coords, aux, coeff);

		BaseFunct[0] = BaseFunct_ID;
		((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta, LineQuadFormula, Needs2ndDer);

		origvaluesD0 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
		origvaluesD1 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);
		origvaluesD2 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D2);

		memset(LocMatrixA, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocMatrixK, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocRhs, 0, N_BaseFunct * SizeOfDouble);

		DOF = GlobalNumbers + BeginIndex[i];
		Vetrex1 = Cell->GetVertex(0);
		Vetrex2 = Cell->GetVertex(1);
		hE = fabs(Vetrex1->GetX() - Vetrex2->GetX());
		//     cout<< " hE  " << hE  <<endl;

		for (j = 0; j < N_Points; j++)
		{
			Coeff = coeff[j];
			c0 = Coeff[0];	   // diffusion
			g0 = Coeff[1];	   // convection
			c1 = Coeff[2];	   // reaction term
			rhsval = Coeff[3]; //rhs

			// double D_L =  TDatabase::ParamDB->REACTOR_P3;
			// double delta0 = TDatabase::ParamDB->DELTA0;
			// double delta1 = TDatabase::ParamDB->DELTA1;

			if (fabs(g0) > 0)
			{

				// Based on Paper
				// https://www.wias-berlin.de/people/john/ELECTRONIC_PAPERS/JKS11.CMAME.pdf

				Pe_K = hE * fabs(g0) / (2. * c0 * order);
				// 	beta based on Lutz book
				//         if(Pe_K>1.)
				//           beta = delta0 * hE/g0;
				//         else
				//           beta = delta1 *hE*hE/D_L ;
				//
				// beta based on 1d Green's formula
				beta = hE * (1. / tanh(Pe_K) - 1. / Pe_K) / (2. * fabs(g0) * order);
			}
			else
			{
				beta = 0.0;
			}

			//  cout<< " Pe_K  " << Pe_K  <<" beta  " << beta  <<endl;
			//  cout<< " c  " << c0 << " g " <<g0 << " c1 " << c1 <<endl;
			Mult = Weights[j] * AbsDetjk[j];
			orgD0 = origvaluesD0[j];
			orgD1 = origvaluesD1[j];
			orgD2 = origvaluesD2[j];

			len += Mult;
			for (k = 0; k < N_BaseFunct; k++)
			{
				test0 = orgD0[k];
				test1 = orgD1[k];

				bgradv = g0 * test1;
				//   	if(k==0)
				// 	   cout<< " bgradv  " << bgradv  <<" bgradv*beta  " << bgradv*beta  <<endl;
				LocRhs[k] += Mult * rhsval * (test0 + beta * bgradv);

				for (l = 0; l < N_BaseFunct; l++)
				{
					ansatz0 = orgD0[l];
					ansatz1 = orgD1[l];
					ansatz2 = orgD2[l];

					val = g0 * ansatz1 * test0; // convective form
					// val  = -g0*ansatz0*test1; // conservative form
					val += c0 * ansatz1 * test1; // difusive term
					val += c1 * ansatz0 * test0;
					//  cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
					LocMatrixA[k * N_BaseFunct + l] += (Mult * val);

					val = (-c0 * ansatz2 + g0 * ansatz1 + c1 * ansatz0) * beta * bgradv;
					LocMatrixK[k * N_BaseFunct + l] += (Mult * val);
				}
			}
		}

		// add to global matrices
		for (j = 0; j < N_BaseFunct; j++)
		{
			TestDOF = DOF[j];
			Rhs[TestDOF] += LocRhs[j];

			begin = RowPtr[TestDOF];
			end = RowPtr[TestDOF + 1];
			for (k = begin; k < end; k++)
			{
				for (l = 0; l < N_BaseFunct; l++)
				{
					if (KCol[k] == DOF[l])
					{
						ValuesA[k] += LocMatrixA[j * N_BaseFunct + l];
						ValuesK[k] += LocMatrixK[j * N_BaseFunct + l];
						break;
					}
				} // for(m=0;m<N_BaseFunct_low
			}	  // for(n=begin;n<end;n++)
		}		  // for(l=0;l<N_BaseFunct_low
	}			  // for(i=0; i<N_Cells_Internal; i++)

	// cout<< " len " << len << endl;
	// exit(0);

	//  // print matrix
	//    for(j=0;j<N_Dof;j++)
	//     {
	//      begin = RowPtr[j];
	//      end = RowPtr[j+1];
	//      for(k=begin;k<end;k++)
	//       {
	//  //       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
	//        cout << "K(" << j << ", "<< KCol[k] << ") = " << ValuesK[k] <<endl;
	//       }

	//  //       cout << "f: " << rhs[j ] <<endl;
	//      cout<<endl;
	//     }
	//  exit(0);
	for (i = 0; i < MaxN_QuadPoints_1D; i++)
	{
		delete[] aux[i];
		delete[] coeff[i];
	}

	delete[] aux;
	delete[] coeff;

	for (i = 0; i < N_Coord + 1; i++)
		delete[] Coords[i];

	delete[] Coords;

} // void TSystemADI1D::AssembleARhs(d

// /**  Assembling A matrix */
void TSystemCD1D::AssembleARhs_DG()
{
	int i, j, k, l, N_Cells_Internal, N_BaseFunct;
	int N_Points, N_Sets = 1, *GlobalNumbers, *BeginIndex, *DOF, *NeibDOF;
	int TestDOF, begin, end, *RowPtr, *KCol;
	int m, N_Joints, Neigh_i;

	double rec_detjk, Neigh_rec_detjk, Neigh_N_BaseFunct;
	double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
	double LocMatrixA11[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocMatrixA12[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocMatrixA21[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocMatrixA22[MaxN_BaseFunctions1D * MaxN_BaseFunctions1D];
	double LocRhs[MaxN_BaseFunctions1D];
	double BdValues[3];
	double **origvaluesD0, **origvaluesD1, Mult;
	double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
	double c0, c1, g0, rhsval, val, len = 0.;
	double **aux, **coeff, *Coeff, **Coords;
	double JointValues[MaxN_BaseFunctions1D], Neigh_JointValues[MaxN_BaseFunctions1D];
	double JointValuesX[MaxN_BaseFunctions1D], Neigh_JointValuesX[MaxN_BaseFunctions1D];
	double Epsilon, Sigma0, Sigma1, h_max, h1, h2;

	bool Needs2ndDer[1], UpdateEdge = FALSE;

	Epsilon = TDatabase::ParamDB->DG_P0;
	Sigma0 = TDatabase::ParamDB->DG_P1;
	Sigma1 = TDatabase::ParamDB->DG_P2;

	TCollection *Coll;
	TBaseCell *Cell, *Neigh = NULL;
	FE1D FEId, Neigh_FEId;
	TFE1D *Element, *Neigh_Element;
	TBaseFunct1D *bf, *Neigh_bf;
	QuadFormula1D LineQuadFormula;
	TQuadFormula1D *qf1;
	TRefTrans1D *F_K;
	BaseFunct1D BaseFunct_ID, BaseFunct[1];

	Coll = FESpace1D->GetCollection();
	GlobalNumbers = FESpace1D->GetGlobalNumbers();
	BeginIndex = FESpace1D->GetBeginIndex();

	RowPtr = A_Intl->GetRowPtr();
	KCol = A_Intl->GetKCol();
	ValuesA = A_Intl->GetEntries();

	N_Cells_Internal = Coll->GetN_Cells();
	Needs2ndDer[0] = FALSE;

	aux = new double *[MaxN_QuadPoints_1D];
	coeff = new double *[MaxN_QuadPoints_1D];

	Coords = new double *[N_Coord + 1];

	for (i = 0; i < N_Coord + 1; i++)
		Coords[i] = new double[MaxN_QuadPoints_1D];

	for (i = 0; i < MaxN_QuadPoints_1D; i++)
	{
		for (int ii = 0; ii < N_Coord; ++ii)
			Coords[ii][i] = 0;

		aux[i] = new double[6];
		coeff[i] = new double[6];
	}

	Intl_L = Coords[N_Coord];

	double BdVal[2];
	BundValues(0, NULL, BdVal);

	// associate each cell with her number in the collection
	for (i = 0; i < N_Cells_Internal; i++)
	{
		Cell = Coll->GetCell(i);
		Cell->SetClipBoard(i);
	}

	for (i = 0; i < N_Cells_Internal; i++)
	{
		UpdateEdge = FALSE;
		Cell = Coll->GetCell(i);
		FEId = FESpace1D->GetFE1D(i, Cell);
		Element = TFEDatabase2D::GetFE1D(FEId);
		bf = Element->GetBaseFunct1D();
		N_BaseFunct = Element->GetN_DOF();
		BaseFunct_ID = Element->GetBaseFunct1D_ID();
		h1 = Cell->GetMeasure();

		//     cout << " BaseFunct_ID " << BaseFunct_ID << endl;
		l = bf->GetPolynomialDegree();
		LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(20 * l);
		qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
		qf1->GetFormulaData(N_Points, Weights, zeta);

		F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
		((TLineAffin *)F_K)->SetCell(Cell);
		((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);
		rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

		//     Bilinear(N_Points, Space_X, Space_Y, Intl_L, aux, coeff);
		// Coords[0][0] = i;
		Bilinear(N_Points, N_Coord + 1, Coords, aux, coeff);

		BaseFunct[0] = BaseFunct_ID;
		((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta, LineQuadFormula, Needs2ndDer);

		origvaluesD0 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
		origvaluesD1 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

		memset(LocMatrixA11, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocMatrixA12, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocMatrixA21, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocMatrixA22, 0, N_BaseFunct * N_BaseFunct * SizeOfDouble);
		memset(LocRhs, 0, N_BaseFunct * SizeOfDouble);

		DOF = GlobalNumbers + BeginIndex[i];

		for (j = 0; j < N_Points; j++)
		{
			Coeff = coeff[j];
			c0 = Coeff[0];	   // diffusion
			g0 = Coeff[1];	   // convection in z direction
			c1 = Coeff[2];	   // reaction term
			rhsval = Coeff[3]; //rhs

			//  cout<< " diffusion  " << c0 << " g0 " <<g0 << " rhsval " << rhsval <<endl;
			Mult = Weights[j] * AbsDetjk[j];
			orgD0 = origvaluesD0[j];
			orgD1 = origvaluesD1[j];

			len += Mult;
			for (k = 0; k < N_BaseFunct; k++)
			{
				test0 = orgD0[k];
				test1 = orgD1[k];

				LocRhs[k] += Mult * rhsval * test0;

				for (l = 0; l < N_BaseFunct; l++)
				{
					ansatz0 = orgD0[l];
					ansatz1 = orgD1[l];

					val = g0 * ansatz1 * test0; // dg non convective form
					// val  = -g0*ansatz0*test1; // dg convective form
					val += c0 * ansatz1 * test1; // difusive term
					val += c1 * ansatz0 * test0;
					// if(i==0)
					// cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
					LocMatrixA11[k * N_BaseFunct + l] += (Mult * val);
				}
			}
		}

		N_Joints = Cell->GetN_Edges();
		for (m = 0; m < N_Joints; m++)
		{
			Neigh = (Cell->GetJoint(m))->GetNeighbour(Cell);

			if (m == 0 && Neigh) // inner left joint, since in 1D we have only 2 joints (really vertices)
			{
				UpdateEdge = TRUE;
				// only first joint (really vertices) will be updated,
				// other joint will be the first joint of the next cell
				h2 = Neigh->GetMeasure();
				h_max = MAX(h1, h2);

				// Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
				// Bilinear(1, N_Coord+1, Coords, aux, coeff);

				// c0 = coeff[0][0]; // no diffusion

				//find the current cell basis value at this joint (X_n+)
				bf->GetDerivatives(D0, -1., JointValues);
				// bf->GetDerivatives(D1, -1., JointValuesX);

				// for(k=0;k<N_BaseFunct;k++)
				//  JointValuesX[k] *=rec_detjk;

				//find the neib cell basis value at this joint (X_n-)
				Neigh_i = Neigh->GetClipBoard();
				Neigh_FEId = FESpace1D->GetFE1D(Neigh_i, Neigh);
				Neigh_Element = TFEDatabase2D::GetFE1D(Neigh_FEId);
				Neigh_bf = Neigh_Element->GetBaseFunct1D();
				Neigh_N_BaseFunct = Neigh_Element->GetN_DOF();

				NeibDOF = GlobalNumbers + BeginIndex[Neigh_i];

				F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
				((TLineAffin *)F_K)->SetCell(Neigh);
				Neigh_rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

				Neigh_bf->GetDerivatives(D0, 1., Neigh_JointValues);
				// Neigh_bf->GetDerivatives(D1, 1., Neigh_JointValuesX);

				// for(k=0;k<Neigh_N_BaseFunct;k++)
				//  Neigh_JointValuesX[k] *=Neigh_rec_detjk;

				//(X_n+, X_n+) (test, ansatz)
				for (k = 0; k < N_BaseFunct; k++)
					for (l = 0; l < N_BaseFunct; l++)
					{
						val = -0.5 * g0 * JointValues[l] * JointValuesX[k] - 0.5 * c0 * Epsilon * JointValues[l] * JointValuesX[k] + (Sigma0 / h_max) * JointValues[l] * JointValues[k] + (Sigma1 / h_max) * JointValuesX[l] * JointValuesX[k];

						//dg convective form
						// val = g0*JointValues[l]*JointValues[k];

						LocMatrixA11[k * N_BaseFunct + l] += val;
						//          cout<<  " val  " << h_max*val  <<endl;
					}

				// e_n
				for (k = 0; k < N_BaseFunct; k++)			// own test X_n+
					for (l = 0; l < Neigh_N_BaseFunct; l++) // neib ansatz X_n-
					{
						val = 0.5 * c0 * Neigh_JointValuesX[l] * JointValues[k] + 0.5 * c0 * Epsilon * Neigh_JointValues[l] * JointValuesX[k] - (Sigma0 / h_max) * Neigh_JointValues[l] * JointValues[k] - (Sigma1 / h_max) * Neigh_JointValuesX[l] * JointValuesX[k];

						//dg convective form
						// val += 0.5*g0*Neigh_JointValues[l]*JointValues[k];

						LocMatrixA21[k * N_BaseFunct + l] += val;
						//            cout<<  " val  " << h_max*val  <<endl;
					}

				// d_n
				for (k = 0; k < Neigh_N_BaseFunct; k++) // neib test X_n-
					for (l = 0; l < N_BaseFunct; l++)	// own ansatz X_n+
					{
						val = -0.5 * c0 * JointValuesX[l] * Neigh_JointValues[k] - 0.5 * c0 * Epsilon * JointValues[l] * Neigh_JointValuesX[k] - (Sigma0 / h_max) * JointValues[l] * Neigh_JointValues[k] - (Sigma1 / h_max) * JointValuesX[l] * Neigh_JointValuesX[k];
						//dg convective form
						// val -= 0.5*g0*JointValues[l]*Neigh_JointValues[k];

						LocMatrixA12[k * N_BaseFunct + l] += val;
						//            cout<<  " val  " << h_max*val  <<endl;
					}

				//(X_n-, X_n-) (test, ansatz)
				for (k = 0; k < Neigh_N_BaseFunct; k++)
					for (l = 0; l < Neigh_N_BaseFunct; l++)
					{
						val = -0.5 * c0 * Neigh_JointValuesX[l] * Neigh_JointValues[k] + 0.5 * c0 * Epsilon * Neigh_JointValues[l] * Neigh_JointValuesX[k] + (Sigma0 / h_max) * Neigh_JointValues[l] * Neigh_JointValues[k] + (Sigma1 / h_max) * Neigh_JointValuesX[l] * Neigh_JointValuesX[k];
						//dg convective form
						// val -= 0.5*g0*Neigh_JointValues[l]*Neigh_JointValues[k];

						LocMatrixA22[k * N_BaseFunct + l] += val;
						//          cout<<  " val  " << h_max*val  <<endl;
					}
			}				 //if(neigh)
			else if (!Neigh) // boundary joint
			{
				if (m == 0) // L_min
				{

					if (cond_Lmin == DIRICHLET)
					{
						//find the current cell basis value at this joint (X_n+)
						bf->GetDerivatives(D0, -1., JointValues);
						bf->GetDerivatives(D1, -1., JointValuesX);

						//  Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
						//            Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
						Bilinear(1, N_Coord + 1, Coords, aux, coeff);

						c0 = coeff[0][0]; // diffusion

						//            BDValue_LMin(Intl_L[0], 0, BdValues);

						BdValues[0] = BdVal[0];

						for (k = 0; k < N_BaseFunct; k++)
							JointValuesX[k] *= rec_detjk;

						for (k = 0; k < N_BaseFunct; k++)
						{
							val = -Epsilon * c0 * JointValuesX[k] * BdValues[0];
							val += (Sigma0 / h1) * JointValues[k] * BdValues[0];

							LocRhs[k] += val;
							// cout<<  " val    rhsval " << val  <<endl;
							for (l = 0; l < N_BaseFunct; l++)
							{
								val = c0 * JointValuesX[l] * JointValues[k] - c0 * Epsilon * JointValues[l] * JointValuesX[k] + (Sigma0 / h1) * JointValues[l] * JointValues[k] + (Sigma1 / h1) * JointValuesX[l] * JointValuesX[k];

								LocMatrixA11[k * N_BaseFunct + l] += val;
								//                 cout<<  " val L_min  " << val  <<endl;
							}
						}
					} //if(cond_Lmin==DIRICHLET)
				}	  // if(m==0)
				else  // L_Max
				{
					// if(cond_Lmax==DIRICHLET || cond_Lmax==NEUMANN)
					if (cond_Lmax == DIRICHLET)
					{
						//find the current cell basis value at this joint (X_n+)
						bf->GetDerivatives(D0, 1., JointValues);
						bf->GetDerivatives(D1, 1., JointValuesX);

						// Intl_L[0] = TDatabase::ParamDB->REACTOR_P13;
						//             Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
						Bilinear(1, N_Coord + 1, Coords, aux, coeff);

						c0 = coeff[0][0]; // diffusion

						// BDVal_LMax(Intl_L[0], 0, BdValues);
						BdValues[0] = BdVal[1];
						BdValues[1] = BdVal[1];

						for (k = 0; k < N_BaseFunct; k++)
							JointValuesX[k] *= rec_detjk;

						for (k = 0; k < N_BaseFunct; k++)
						{
							val = Epsilon * c0 * JointValuesX[k] * BdValues[0];
							val += (Sigma1 / h1) * JointValuesX[k] * BdValues[1];

							LocRhs[k] += val;
							//                 cout<<  " val    rhsval " << val  <<endl;
							for (l = 0; l < N_BaseFunct; l++)
							{
								val = -c0 * JointValuesX[l] * JointValues[k] + c0 * Epsilon * JointValues[l] * JointValuesX[k] + (Sigma0 / h1) * JointValues[l] * JointValues[k] + (Sigma1 / h1) * JointValuesX[l] * JointValuesX[k];

								LocMatrixA11[k * N_BaseFunct + l] += val;
								//                cout<<  " val L_max " << val  <<endl;
							}
						}
					}
				}
			}
		}

		// add to global matrices
		for (j = 0; j < N_BaseFunct; j++)
		{
			TestDOF = DOF[j];
			Rhs[TestDOF] += LocRhs[j];

			begin = RowPtr[TestDOF];
			end = RowPtr[TestDOF + 1];
			for (k = begin; k < end; k++)
			{
				for (l = 0; l < N_BaseFunct; l++)
				{
					if (KCol[k] == DOF[l])
					{
						ValuesA[k] += LocMatrixA11[j * N_BaseFunct + l];
						//             cout << TestDOF << " " << KCol[k] << endl;
						break;
					}
				} // for(m=0;m<N_BaseFunct

				// add edge integrals
				if (UpdateEdge)
				{
					for (l = 0; l < Neigh_N_BaseFunct; l++)
					{
						if (KCol[k] == NeibDOF[l])
						{
							//             cout << TestDOF << " " << KCol[k] << endl;
							ValuesA[k] += LocMatrixA21[j * N_BaseFunct + l];
							break;
						}
					} // for(l=0;l<Neigh_N_BaseFunct;l++)
				}

			} // for(n=begin;n<end;n++)
		}	  //  for(j=0;j<N_BaseFunct;j++)

		// add edge integrals
		if (UpdateEdge)
		{
			for (j = 0; j < Neigh_N_BaseFunct; j++)
			{
				TestDOF = NeibDOF[j];

				begin = RowPtr[TestDOF];
				end = RowPtr[TestDOF + 1];
				for (k = begin; k < end; k++)
				{
					for (l = 0; l < N_BaseFunct; l++)
					{
						if (KCol[k] == DOF[l])
						{
							ValuesA[k] += LocMatrixA12[j * N_BaseFunct + l];
							//             cout << TestDOF << " " << KCol[k] << endl;
							break;
						}
					} //  for(l=0;l<N_BaseFunct;l++)

					for (l = 0; l < Neigh_N_BaseFunct; l++)
					{
						if (KCol[k] == NeibDOF[l])
						{
							ValuesA[k] += LocMatrixA22[j * N_BaseFunct + l];
							//             cout << TestDOF << " " << KCol[k] <<  " " <<  LocMatrixA22[j*N_BaseFunct + l] << endl;
							break;
						}
					} //  for(l=0;l<N_BaseFunct;l++)
				}	  // for(k=begin;k<end;k++)
			}		  // for(j=0;j<Neigh_N_BaseFunct;j++)
		}			  // if(UpdateEdge)

	} // for(i=0; i<N_Cells_Internal; i++)

	// cout<< " len " << len << endl;
	//  exit(0);
	//  print matrix
	//     for(j=0;j<N_Dof;j++)
	//      {
	//  //      begin = RowPtr[j];
	//  //      end = RowPtr[j+1];
	//  //      for(k=begin;k<end;k++)
	//  //       {
	//  //        cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
	//  //       }

	//         cout << "f: " << rhs[j ] <<endl;
	//       cout<<endl;
	//       }
	//   exit(0);
	for (i = 0; i < MaxN_QuadPoints_1D; i++)
	{
		delete[] aux[i];
		delete[] coeff[i];
	}

	delete[] aux;
	delete[] coeff;
} // void TSystemCD1D::AssembleARhs(d

// Print the Solution to the terminal Output
void TSystemCD1D::printSolution()
{
	int N_DOF = this->GetN_Dof();
	double *Solution = this->Sol;
	for (int i = 0; i < N_DOF; i++)
		cout << Solution[i] << ",";
	cout << endl;
}

void TSystemCD1D::generateVTK()
{
	std::string filename = TDatabase::ParamDB->VTKBASENAME;
	int N_DOF = this->GetN_Dof();
	double start = TDatabase::ParamDB->START_X;
	double end = TDatabase::ParamDB->END_X;
	int N_Cells = this->GetN_Cells();
	double h = (end - start) / N_Cells;
	int N_Nodes_Element = TDatabase::ParamDB->ANSATZ_ORDER + 1;
	double *Solution = this->Sol;

	filename.append(".vtk");
	std::ofstream myfile(filename);

	myfile.flags(std::ios::dec | std::ios::scientific);
	myfile.precision(6);

	// LIne 1 - VTK Version
	myfile << "# vtk DataFile Version 1.0" << std::endl;

	//Line 2 - Title
	myfile << " Solution 1D " << std::endl
		   << std::endl;

	//Line 3 - Data type ( ASCII / BINARY )
	myfile << "ASCII" << std::endl;

	//Line 3 - Structured or unstructured grid ( STRUCTURED_GRID / UNSTRUCTURED_GRID )
	myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

	//Line 4 - no of data points :: syntax -> POINTS <no of points> <datatype>
	myfile << "POINTS " << N_DOF << " float" << std::endl;

	for (int i = 0; i < N_DOF; i++)
		myfile << i * (h) << " "
			   << "0.0000000 "
			   << "0.000000" << std::endl;

	myfile << std::endl;

	// CELLS -> syntax : CELLS <no of cells> <no of parameters  totally needed to define the total cell points>
	// for eg: in 1D total points in cells is 2 , So the last parameter will be (2+1)* No of cells

	myfile << "CELLS " << N_Cells << " " << (N_Nodes_Element + 1) * N_Cells << std::endl;

	int node = 0;
	for (int i = 0; i < N_Cells; i++)
		myfile << "2 " << node << " " << ++node << std::endl;

	myfile << std::endl;

	// CELL TYPES : syntax: CELL_TYPES <No of cells>
	// cell type for all the cells,
	// for 1d Element - cell type is 3

	myfile << "CELL_TYPES " << N_Cells << std::endl;

	for (int i = 0; i < N_Cells; i++)
		myfile << "3 ";

	myfile << std::endl
		   << std::endl
		   << std::endl;

	// POINT DATA : syntax - PONT_DATA <no of points>
	// < scalar or vector> < Datatype>
	//"LOOKUPTABLE" < lookuptable type >

	myfile << "POINT_DATA " << N_DOF << std::endl;
	myfile << "SCALARS "
		   << "1D_Solution "
		   << "float" << std::endl;
	myfile << "LOOKUP_TABLE "
		   << "default" << std::endl;

	for (int i = 0; i < N_DOF; i++)
		myfile << Solution[i] << std::endl;

	myfile.close();

	std::cout << "VTK File " << filename << " has been generated" << std::endl;
}

void TSystemCD1D::plotGNU()
{

	int N_DOF = this->GetN_Dof();
	double start = TDatabase::ParamDB->START_X;
	double end = TDatabase::ParamDB->END_X;
	int N_Cells = this->GetN_Cells();
	double h = (end - start) / N_Cells;
	double *Solution = this->Sol;
	std::ofstream file;

	file.open("GNUPLOT_FILE.dat");

	//file.flags( std::ios::dec | std::ios::scientific);
	file.precision(6);

	file << "#"
		 << "   "
		 << "X"
		 << "   "
		 << "Y" << std::endl;

	for (int i = 0; i < N_DOF; i++)
	{
		file << i * h << "   " << Solution[i] << std::endl;
	}

	file.close();

	system("gnuplot -e \"plot 'GNUPLOT_FILE.dat' using 1:2 with linespoints; pause -1 ; exit\" ");
}

TSystemCD1D::~TSystemCD1D()
{
	//  delete [] Nodal_sol;
	//  delete [] oldsol;
	//  delete [] rhs;
	//  delete [] sol;
	//  delete [] B;
	//  delete [] defect;
}
