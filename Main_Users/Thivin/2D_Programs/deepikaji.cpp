#include <Domain.h>
#include <Mapper.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <LinAlg.h> 
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <NSE2D_ParamRout.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <SystemNSE2D.h>
#include <Assemble2D.h>




#include <vector>
#include<fstream>

using namespace std;

void Assembly_poisson_2D(double quad_wt, double *coeff, double *param,
                                      double hK, double **derivatives, int *N_BaseFuncts,double ***LocMatrices, double **LocRhs)
{

  double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2];
  double **K11, **K12, **K21, **K22, *F1, *F2;
  K11 = LocMatrices[0];
  K12 = LocMatrices[1];
  K21 = LocMatrices[2];
  K22 = LocMatrices[3];

  F1 = LocRhs[0];
  F2 = LocRhs[1];
  double mu = 1;

  for (int i = 0; i < N_BaseFuncts[0]; i++){
    for (int j = 0; j < N_BaseFuncts[0]; j++){
     

       K11[i][j] += quad_wt * (N[i]*N[j]); 

      //cout << "Nx[ " << i << "] = " << Nx[i] << endl;

      K12[i][j] += 0.0;
      K21[j][i] += 0.0;
 
      K22[i][j] += quad_wt * (N[i]*N[j]); 

      /* NON LINEAR PART */
    }
    /* RHS */

    F1[i] = 0.;

    F2[i] = 0.;
  }

}



void BilinearCoeffs(int n_points, double *x, double *y,double **parameters, double **coeffs)
{
  double angle = 0, v1, v2;
  int i;
  double *coeff;

  //v1 = cos(angle);
  //v2 = sin(angle);
}


void BoundCondition_Y(int i, double t, BoundCond &cond)
{
      cond = DIRICHLET;

}

void BoundCondition_X(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

double boundval = 0;
void BoundValue_X(int BdComp, double Param, double &value)
{
    value = 0.;
}


void BoundValue_Y(int BdComp, double Param, double &value)
{
  //cout << "Bd Comp = " << BdComp << ", Param = " << Param << ", value = " << value << endl;
	value = 0;
}

typedef void AssembleFct2D(double, double *, double, double **, 
                           int *, double ***, double **);


void NSType1Galerkin_deepikaJi(double Mult, double *coeff,
                    double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test00*ansatz00);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// void NSType1Upwind (double Mult, double *coeff,
//                     double *param, double hK,
//                     double **OrigValues, int *N_BaseFuncts,
//                     double ***LocMatrices, double **LocRhs)
// {
//   double **MatrixA, **MatrixB1, **MatrixB2;
//   double *Rhs1, *Rhs2, val;
//   double *MatrixRow, *MatrixRow1, *MatrixRow2;
//   double ansatz10, ansatz01;
//   double test00, test10, test01;
//   double *Orig0, *Orig1, *Orig2, *Orig3;
//   int i,j, N_U, N_P;
//   double c0, c1, c2;

//   MatrixA = LocMatrices[0];
//   MatrixB1 = LocMatrices[1];
//   MatrixB2 = LocMatrices[2];

//   Rhs1 = LocRhs[0];
//   Rhs2 = LocRhs[1];

//   N_U = N_BaseFuncts[0];
//   N_P = N_BaseFuncts[1];

//   Orig0 = OrigValues[0];         // u_x
//   Orig1 = OrigValues[1];         // u_y
//   Orig2 = OrigValues[2];         // u
//   Orig3 = OrigValues[3];         // p

//   c0 = coeff[0];                 // nu
//   c1 = coeff[1];                 // f1
//   c2 = coeff[2];                 // f2

//   //  cout << "c3";

//   for(i=0;i<N_U;i++)
//   {
//     MatrixRow = MatrixA[i];
//     test10 = Orig0[i];
//     test01 = Orig1[i];
//     test00 = Orig2[i];

//     Rhs1[i] += Mult*test00*c1;
//     Rhs2[i] += Mult*test00*c2;

//     for(j=0;j<N_U;j++)
//     {
//       ansatz10 = Orig0[j];
//       ansatz01 = Orig1[j];

//       val = c0*(test10*ansatz10+test01*ansatz01);
//       //    val += c3 * Orig2[j] *test00;

//       MatrixRow[j] += Mult * val;
//     }                            // endfor j
//   }                              // endfor i

//   for(i=0;i<N_P;i++)
//   {
//     MatrixRow1 = MatrixB1[i];
//     MatrixRow2 = MatrixB2[i];

//     test00 = Orig3[i];

//     for(j=0;j<N_U;j++)
//     {
//       ansatz10 = Orig0[j];
//       ansatz01 = Orig1[j];

//       val = -Mult*test00*ansatz10;
//       MatrixRow1[j] += val;

//       val = -Mult*test00*ansatz01;
//       MatrixRow2[j] += val;
//     }                            // endfor j

//   }                              // endfor i
// }



int main (int argc, char** argv)
{

    int N_Cells, ORDER, N_U, N_P, N_DOF, N_Active;
	double *sol, *sol_buffer, *rhs, *defect, t1, t2, errors[4], residual, impuls_residual;
	double limit, u_error[4], p_error[2];
	srand(time(NULL));
	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();

    /* ---------------  Collection of all the cells --------------/
	* -- Connectivity Information , Global to local Node NUmbering -- */
	TCollection *coll;
    TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2], *fesprhs[3];
    TFEVectFunct2D *Velocity;
    TFEFunction2D *u1, *u2, *Pressure, *fefct[3];
    TOutput2D *Output;
    TSystemNSE2D *SystemMatrix;


    const char vtkdir[] = "VTK";
    char *PsBaseName, *VtkBaseName, *GEO;
    char UString[] = "u";
    char PString[] = "p";

    std::ostringstream os;
    os << " ";


    Domain = new TDomain(argv[1]);  

    OpenFiles();
    OutFile.setf(std::ios::scientific);

    Database->CheckParameterConsistencyNSE();
    Database->WriteParamDB(argv[0]);



    if(TDatabase::ParamDB->MESH_TYPE==0)
   {
    GEO = TDatabase::ParamDB->GEOFILE;
    Domain->Init(NULL, GEO);
   }
     else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {
         cout << " MESH - GMSH GEN " << endl;
       Domain->GmshGen(TDatabase::ParamDB->GEOFILE); 
       OutPut("GMSH used for meshing !!!" << endl);
    }
       else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     } 

    if(TDatabase::ParamDB->MESH_TYPE==0){
        for(int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
            Domain->RegRefineAll();  
    }

    if(TDatabase::ParamDB->WRITE_PS)
    {
    // write grid into an Postscript file
        os.seekp(std::ios::beg);
        os << "Domain" << ".ps" << ends;
        Domain->PS(os.str().c_str(),It_Finest,0);
    }

    if(TDatabase::ParamDB->WRITE_VTK)
        { mkdir(vtkdir, 0777); }  


    int VELOCITY_ORDER  = 1;
    int PRESSURE_ORDER = 1;

    coll=Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells : " << N_Cells <<endl);

    char nameString[] = "name";
    char uString[] = "u";
    char pString[] = "p";


    Velocity_FeSpace = new TFESpace2D(coll,nameString,uString,BoundCondition_X,VELOCITY_ORDER, NULL);

    Pressure_FeSpace = new TFESpace2D(coll,nameString,pString,BoundCondition_X,PRESSURE_ORDER, NULL);
    int pressure_space_code = PRESSURE_ORDER;
    int velocity_space_code = VELOCITY_ORDER;
    N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
    N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();    
    int N_TotalDOF = 2*N_U + N_P;

    OutPut("Dof Velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("Dof Pressure : "<< setw(10) << N_P << endl);
    OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);


    //======================================================================
    // construct all finite element functions
    //======================================================================
    sol = new double[N_TotalDOF];
    rhs = new double[N_TotalDOF];

    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
    Pressure = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);


    fesp[0] = Velocity_FeSpace;
    fefct[0] = u1;
    fefct[1] = u2;

    TAuxParam2D* auxerror = NULL;

    TAuxParam2D *aux = new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
			   NSN_ParamFctVelo,
                           NSN_FEValuesVelo,
                           fesp, fefct,
                           NSFctVelo,
                           NSFEFctIndexVelo, NSFEMultiIndexVelo,
                           NSN_ParamsVelo, NSBeginParamVelo);


    int NSEType = TDatabase::ParamDB->NSTYPE;


    SystemMatrix = new TSystemNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, GALERKIN, NSEType, DIRECT, NULL,NULL,NULL);

    

    ///////// CODE BLOCK FOR INIT NAVIER STOKES - ///////////////////////////////////

    SystemMatrix->Init(BilinearCoeffs, BoundCondition_X, BoundValue_Y, BoundValue_Y, aux, auxerror);

    TDiscreteForm2D *DiscreteFormGalerkin , *DiscreteFormARhs , *DiscreteFormNL ,*DiscreteFormUpwind;
    
    BoundCondFunct2D *BoundaryConditions[2]; BoundValueFunct2D *BoundaryValues[2];
    CoeffFct2D* LinCoeffs[1];

    // save the boundary condition
    BoundaryConditions[0] = BoundCondition_X; BoundaryConditions[1] = BoundCondition_X;  

    // save the boundary values  
    BoundaryValues[0] = BoundValue_Y;BoundaryValues[1] = BoundValue_Y;
    // save the nse bilinear coefficient   
    LinCoeffs[0] = BilinearCoeffs;


    // DISCRETE FORM PARAMETERS 
    char Galerkin[] = "Galerkin";
    char all[] = "all";
    char UpwindString[] = "Upwind";
    int NSType1N_Terms = 4;
    MultiIndex2D NSType1Derivatives[4] = { D10, D01, D00, D00 };
    int NSType1SpaceNumbers[4] = { 0, 0, 0, 1 };
    int NSType1N_Matrices = 3;
    int NSType1RowSpace[3] =   { 0, 1, 1 };
    int NSType1ColumnSpace[3] = { 0, 0, 0 };
    int NSType1N_Rhs = 2;
    int NSType1RhsSpace[2] = { 0, 0 };
    ManipulateFct2D *manipulate;
    manipulate = NULL;

    // double x;
    DiscreteFormGalerkin = new TDiscreteForm2D(Galerkin, all,
                                                NSType1N_Terms, NSType1Derivatives, 
                                                NSType1SpaceNumbers,
                                                NSType1N_Matrices, NSType1N_Rhs, 
                                                NSType1RowSpace, NSType1ColumnSpace,
                                                NSType1RhsSpace, Assembly_poisson_2D, LinCoeffs, nullptr);



    // // DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, all, 
    // //                                         NSType1N_Terms, NSType1Derivatives, 
    // //                                         NSType1SpaceNumbers,
    // //                                         NSType1N_Matrices, NSType1N_Rhs, 
    // //                                         NSType1RowSpace, NSType1ColumnSpace,
    // //                                         NSType1RhsSpace, NSType1Upwind, LinCoeffs, nullptr);

    // DiscreteFormARhs = DiscreteFormUpwind;     
    // DiscreteFormNL = NULL;

    // // MATRICES:
    // TSquareMatrix2D *SqmatrixA11, *SqmatrixA12, *SqmatrixA21, *SqmatrixA22, *SQMATRICES[9];
    // TMatrix2D *MatrixB1, *MatrixB2, *MatrixB1T, *MatrixB2T, *MATRICES[8];

    // SQMATRICES[0] = SqmatrixA11;
    // MATRICES[0] = MatrixB1;
    // MATRICES[1] = MatrixB2;

    // SQMATRICES[0]->Reset();
    // MATRICES[0]->Reset();
    // MATRICES[1]->Reset();

    // int N_SquareMatrices = 1;
    // int N_RectMatrices = 2;
    // int  N_Rhs = 2;
    // int N_FESpaces = 2;  

    // TFESpace2D *FeSpaces[5];

    // fesprhs[0] = FeSpaces[0];
    // fesprhs[1] = FeSpaces[0];
    // fesprhs[2] = FeSpaces[1];   

    // AMatRhsAssemble = new TAssembleMat2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES, N_RectMatrices, MATRICES,
    //                           N_Rhs, RHSs, fesprhs, DiscreteFormARhs, BoundaryConditions, BoundaryValues, NSEaux);
    
    // AMatRhsAssemble->Init();  

    //////////////// END OF CODE BLOCK - NAVIER STOKES  - ////////////////////////////////////
    
    

    SystemMatrix->Assemble(sol, rhs);

    // calculate the residual
    defect = new double[N_TotalDOF];
    memset(defect,0,N_TotalDOF*SizeOfDouble);

    SystemMatrix->GetResidual(sol, rhs, defect);

    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    


    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput2D(2, 1, 1, 2, Domain);

    Output->AddFEVectFunct(Velocity);
    Output->AddFEFunction(Pressure);

    int img=1;

    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     } 


    residual =  Ddot(N_TotalDOF, defect, defect);
    impuls_residual = Ddot(2*N_U, defect, defect);  

    OutPut("Nonlinear iteration step   0");
    OutPut(setw(14) << impuls_residual);
    OutPut(setw(14) << residual-impuls_residual);
    OutPut(setw(14) << sqrt(residual) << endl);  

    SystemMatrix->Solve(sol, rhs);


    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
      residual =  Ddot(N_TotalDOF, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect); 

    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20FEFunction(sol+2*N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

    
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput2D(2, 1, 1, 2, Domain);

    Output->AddFEVectFunct(Velocity);
    Output->AddFEFunction(Pressure);
   
  //  u1->Interpolate(ExactU1);
  //  u2->Interpolate(ExactU2);
  //  Pressure->Interpolate(ExactP);
   
    if(TDatabase::ParamDB->WRITE_VTK)
        {
        os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
            else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
            else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
        }   


    CloseFiles();



    return 0;
}   