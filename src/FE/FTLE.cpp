#include<AllClasses.h>
#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <DeformMesh2D.h>
#include <LinAlg.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include<FEVectFunct2D.h>
#include <Constants.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <NodalFunctional2D.h>
#include<BaseCell.h>
#include<Collection.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include<vector>

#include<FTLE.h>


// COnstructor
FTLE::FTLE(TFEVectFunct2D* FEVectFunction,int searchDepth)
{   
    TFESpace2D* fespace =  FEVectFunction->GetFESpace2D();
    TCollection* coll =  fespace->GetCollection();
    VelocityVectFunction = FEVectFunction;
    int N_DOF = fespace->GetN_DegreesOfFreedom();
    int N_cells =  fespace->GetN_Cells();
    int N_Points;
    double *xi, *eta;
    RefTrans2D RefTrans, *RefTransArray;

    // Construct the arrays for position and the cellId's
    x_pos =  new double[N_DOF]();
    y_pos = new double[N_DOF]();

    // Updates the position of All particles in the domain. 
    fespace->GetDOFPosition(x_pos,y_pos,cellId);

    int* GlobalDOFArray = fespace->GetGlobalNumbers();
    int* BeginIndexArray = fespace->GetBeginIndex();

    cout << " N_cells : " << N_cells <<endl;
    cout << " N_DOF : " << N_DOF <<endl;    
    //Create a vector of vector with size N_DOF;
    
   N_DOF = 100;
    std::vector< std::vector<int> > mapping(100,std::vector<int>());

    for( int k = 0 ; k < N_DOF; k++)
    {
        mapping[k].push_back(0);
    }
    
    exit(0);


    for ( int cellno = 0 ; cellno < N_cells ; cellno++)
    {
        int* globalDOFCell = GlobalDOFArray + BeginIndexArray[cellno];
        // Get the cell , Reference transformation and its basis functions to get the local DOF> 
        TBaseCell *cell = coll->GetCell(cellno);
        FE2D FEid = fespace->GetFE2D(cellno, cell); 
        RefTrans = RefTransArray[FEid];
        BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEid);

        TNodalFunctional2D* nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(FEid);
        nf->GetPointsForAll(N_Points, xi, eta);
        // printf(" hello");
        for ( int j = 0  ; j < N_Points ; j++)
        {
            int dof = globalDOFCell[j];
            printf("%d ", dof);
            // Cell2DOFMapping[cellno].push_back(dof);
            // printf("%d ", j);
            // DOF TO CELL MAPPING 
            // DOF2CellMapping[dof].push_back(cellno);
            // printf("%d ", j);
        }
        
        
    }

}

 double* FTLE::obtainVelocity()
 {
     double* uVector = VelocityVectFunction->GetComponent(0)->GetValues();
     double* vVector = VelocityVectFunction->GetComponent(1)->GetValues();


    for (int i = 0 ; i < N_Particles; i++)
    {
        int currentCell = cellId[i];
        
        //get all the Cells that are neibhour to those current Vertex. 


    }


 }