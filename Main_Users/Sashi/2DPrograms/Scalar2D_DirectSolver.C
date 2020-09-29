// =======================================================================
//
// Purpose:     main program with parallel solver (no multigrid solver)
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.07.2009
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>

double bound = 0;

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <MainUtilities.h>
#include <Upwind.h>
#include <FixedPointIte.h>

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

// ======================================================================
// include the required example file
// ======================================================================
// #include "../Examples/ConvDiff2D/ReactionDominate.h" 
// #include "../Examples/ConvDiff2D/Plane.h"  // unit square
// #include "../Examples/ConvDiff2D/Hemker.h" // circle in a channel
#include "../Examples/CD_2D/TwoBoundaryLayers.h"  // unit square
// #include "../Examples/ConvDiff2D/Smooth.h"  // unit square
// #include "../Examples/CD_2D/TwoInteriorLayers.h"  // unit square
// #include "../Examples/ConvDiff2D/SineLaplaceDiriHom.h" 
// #include "../Examples/ConvDiff2D/SineLaplace.h" 

// ======================================================================
// subroutines for the main program
// ======================================================================
// #include "../2DSubPrograms/Scalar2D_par_sub.h" 


// ======================================================================
// utilities for main program
// ======================================================================
void SetCellDataForDepBasis(int N_Cells, TCollection *coll)
{
 int i;
 
 double Q1, Q2;
 
 TBaseCell *cell;   
   
  for(i=0; i<N_Cells; i++)
   {
    cell =  coll->GetCell(i);
    
 
    GetQ1Q2(Q1, Q2);
     
    cell->SetDepBasisValues();
     
     
   }
   
   
 }




// ======================================================================
// main program
// ======================================================================
int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll;
  TBaseCell *cell;
  TFESpace2D *scalar_space;
  TOutput2D *Output;
  TAuxParam2D *aux;
  TSquareStructure2D *sqstructureA;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[1];
  TSquareMatrix2D **MatricesA, *Tmp_MatricesA;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;

  TDiscreteForm2D *DiscreteForm;
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormUpwind;
  TFEFunction2D *C;
  TFESpace2D *fesp[2], *ferhs[3];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };

  int i,j,k,l,ret, ORDER, N_DOF, N_Active, N_NonActive;
  int N_Unknowns, img=1;
  int N_LinIter, LEVELS, N_Cells;

  double solver_time, hmin, hmax;
  double *sol, *oldsol, *rhs, *defect, *PostSol;
  double *RHSs[1], t1, t2, errors[4];
  double  *l2, *h1, *sd;

  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  char CString[] = "c";
  char UString[] = "u";
  char PString[] = "p";
  char ReadinDat [] = "readin.dat";
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char Name[] = "Temp";
  char Description[] = "description";
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  char SDFEMString[] = "SDFEM";
  char UpwString[] = "Upwind";

  std::ostringstream os, opts;
  os << " ";
  opts << " ";
//======================================================================
// read parameter file Readin.dat
//======================================================================
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);


  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();

//======================================================================
// copy read parameters into local variables
//======================================================================
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  LEVELS = TDatabase::ParamDB->LEVELS;

  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
//======================================================================
// define discrete form
//======================================================================
   DiscreteFormGalerkin = new TDiscreteForm2D
                 (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace, CD_ColumnSpace, CD_RhsSpace,
                  BilinearAssemble, BilinearCoeffs, NULL);  
  
//======================================================================
// generate mesh and refun (if needed)
//======================================================================
    Domain->Init(PRM, GEO);

//       write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain_Coarse" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);

// 
// /* generate special mesh for Hemker example */
//     if(!strcmp(GEO, "InitGrid"))
//      if(TDatabase::ParamDB->REACTOR_P25)
//          MeshReGen_HemkerResolved(Domain);
//      else
//          MeshReGen_Hemker(Domain);

  // refine grid
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
     Domain->RegRefineAll();


      // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);

//======================================================================
// include the boundary condition and boundary values from the example file
//======================================================================
  BoundCondFunct2D *BoundaryConditions[1] = { BoundCondition };
  BoundValueFunct2D *BoundaryValues[1] = { BoundValue };

//   CoeffFct2D *Coefficients[1];
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

  if(ORDER>=1100)
   TDatabase::ParamDB->DEPENDENT_BASIS = 1;

//======================================================================
// loop over all levels (not a multigrid level but for convergence study)
//======================================================================

  for(i=0;i<LEVELS;i++)
   {
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    solver_time = 0.0;
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    if(i)
     Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    OutPut( "number of cells: " << N_Cells << endl);
    coll->GetHminHmax(&hmin,&hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
    cout << endl << endl;

//======================================================================
// construct all finite element spaces
//======================================================================
    scalar_space =  new TFESpace2D(coll, Name, Description,
                                   BoundCondition, ORDER, NULL);

    N_DOF = scalar_space->GetN_DegreesOfFreedom();
    N_Active = scalar_space->GetActiveBound();
    N_NonActive = N_DOF - N_Active;
    OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
//======================================================================
// construct all finite element functions
//======================================================================
    N_Unknowns = N_DOF;
    sol = new double[N_Unknowns];
    PostSol = new double[N_Unknowns];
    oldsol = new double[N_Unknowns];
    rhs = new double[N_Unknowns];
    defect = new double[N_Unknowns];

    memset(sol, 0, N_Unknowns*SizeOfDouble);
    memset(PostSol, 0, N_Unknowns*SizeOfDouble);
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);
    memset(rhs, 0, N_Unknowns*SizeOfDouble);

    C = new TFEFunction2D(scalar_space, UString, UString, sol, N_DOF);

//======================================================================
// allocate memory for all matrices
//======================================================================
    sqstructureA = new TSquareStructure2D(scalar_space);
    sqstructureA->Sort();

    sqmatrixA = new TSquareMatrix2D(sqstructureA);
//======================================================================
// assemble all matrices
//======================================================================
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);

    fesp[0] = scalar_space;
    ferhs[0] = scalar_space;

    switch(TDatabase::ParamDB->DISCTYPE)
     {
      case GALERKIN:
           DiscreteForm = DiscreteFormGalerkin;
      break;

      case SDFEM:
           DiscreteForm = DiscreteFormSDFEM;
      break;

      case UPWIND:
           DiscreteForm = DiscreteFormUpwind;
      break;

      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }

      // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset();
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    if(TDatabase::ParamDB->DEPENDENT_BASIS)
    {
     SetCellDataForDepBasis(N_Cells, coll);      
    }


    // assemble
    Assemble2D(1, fesp,
               1, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteForm,
               BoundaryConditions,
               BoundaryValues,
               aux);

     delete aux;

// // apply local projection stabilization method
//      if(TDatabase::ParamDB->REACTOR_P19>0)
//       {
//        if(TDatabase::ParamDB->REACTOR_P20==0)
//         UltraLocalProjection(SQMATRICES[0], TRUE);
//        else if(TDatabase::ParamDB->REACTOR_P20==1)
//         StreamlineUltraLocalProjection(SQMATRICES[0], TRUE);
// //        else if(TDatabase::ParamDB->REACTOR_P20==2)
// //         StreamlineUltraLocalProjection(SQMATRICES[0], TRUE);
//       }
// 
//     // set rhs for Dirichlet nodes
//     memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
// 
// //======================================================================
// // solve the system
// //======================================================================
//     // compute defect
//     memset(defect,0,N_Unknowns*SizeOfDouble);
// 
// //     OutPut("norm of solution " <<  sqrt(Ddot(N_Active,sol,sol))  << endl);
//     t1 = GetTime();
//     DirectSolver(sqmatrixA, rhs, sol);
//     t2 = GetTime();
//     OutPut( "time for AMG solving: " << t2-t1 << endl);
// //     OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
// 
// //======================================================================
// // post processing
// //======================================================================
//     AdaptivePostProcess(C, PostSol, TRUE);
// 
// 
// //======================================================================
// // produce outout
// //======================================================================
//     if(TDatabase::ParamDB->WRITE_PS)
//     {
//       // write grid into an Postscript file
//       os.seekp(std::ios::beg);
//       os << PsBaseName << i << ".ps" << ends;
//       Domain->PS(os.str().c_str(),It_Finest,0);
//     }
// 
//     Output = new TOutput2D(2, 2, 1, 1,Domain);
// 
//     Output->AddFEFunction(C);
// 
// //     C->Interpolate(Exact);
// 
//     if(TDatabase::ParamDB->WRITE_VTK)
//      {
// //       os.seekp(std::ios::beg);
// //        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
// //          else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
// //           else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
// //            else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
// //             else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
// //       Output->WriteVtk(os.str().c_str());
// //       img++;
//       Output->ParMooN_WriteVTK(img);
//      }
// 
//     memcpy(sol, PostSol, N_Active*SizeOfDouble);
// 
//     if(TDatabase::ParamDB->WRITE_VTK)
//      {
//       os.seekp(std::ios::beg);
//        if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
//       Output->WriteVtk(os.str().c_str());
// //       img++;
//       Output->ParMooN_WriteVTK(img);
//      }
// 
//     ComputeExtremalValues(N_Unknowns, C,errors);
//     OutPut(setprecision(4) << "C min:= " << errors[0] << ", C max:= " << errors[1] -1<<endl);
// 
//     // measure errors to known solution
//     if(TDatabase::ParamDB->MEASURE_ERRORS)
//      {
//       aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
// 
//       C->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
//                    BilinearCoeffs, aux, 1, fesp, errors);
// 
//       delete aux;
// 
// 
//       l2[i] = errors[0];
//       h1[i] = errors[1];
//       sd[i] = errors[2];
// 
//       if (i)
//        {
//         OutPut( "L2: " << errors[0]<< " order " <<  log(l2[i-1]/l2[i])/ln2 << endl);
//         OutPut( "H1-semi: " << errors[1] << " order " << log(h1[i-1]/h1[i])/ln2 << endl);
//         OutPut( "SD: " << errors[2] << endl);
//        }
//       else
//        {
//         OutPut( "L2: " << errors[0] << endl);
//         OutPut( "H1-semi: " << errors[1] << endl);
//         OutPut( "SD: " << errors[2] << endl);
//        }
// 
// 
//      } // if(TDatabase::ParamDB->MEASURE_ERRORS)
// 

    OutPut("memory after: " << setw(10) << GetMemory() << endl);
   } // for(i=0;i<LEVELS;i++)



  CloseFiles();
  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  return 0;
}
