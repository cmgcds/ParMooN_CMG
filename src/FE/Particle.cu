/** ************************************************************************ 
* @brief     source file for GPU Routines of Particle Deposition within the problem
             
* @brief  Description :

            Contains all the functions that are needed to calculate the particle position within the cells.

            Functions and formulas are same as implemented in the sequential Implementation of ParMooN.
            They are just simplified ( in Datastructure wise ) for implementing in GPU

* @author  Thivin Anandh D
* @date      20-Jun-2023
* @History   Jun 20  - Creation of Cuda Routines
 ************************************************************************  */

#include <cuda.h>
#include "helper_cuda.h"
#include <AllClasses.h>
#include <FEDatabase3D.h>
#include "Particle.h"

void TParticles::SetupCudaDataStructures(TFESpace3D* fespace){
    // get the collection of cells from the fespace
    TCollection *coll = fespace->GetCollection();

    // get No of cells in the collection
    int N_Cells = coll->GetN_Cells();

    // get N_dof
    int N_DOF = fespace->GetN_DegreesOfFreedom();

    // get the global indices and begin index of the cells from fespace
    int *GlobalNumbers = fespace->GetGlobalNumbers();
    int *BeginIndex = fespace->GetBeginIndex();

    // get the last value in the begin index array
    int size_of_global_numbers = BeginIndex[N_Cells];

    // Allocate memory for the cell vertices
    h_m_cell_vertices = new double[4 * N_Cells];

    // Allocate memory for global indices of the cells
    h_m_global_dof_indices = new int[size_of_global_numbers];

    // Allocate memory for the begin index of the cells
    h_m_begin_indices  = new int[N_Cells + 1];

    // allocate memory for velocity values
    h_m_velocityX = new double[N_DOF];
    h_m_velocityY = new double[N_DOF];
    h_m_velocityZ = new double[N_DOF];

    // hardcode the value of n_basis functions based on the first cell
    // Get the first cell
    TBaseCell *cell_0 = coll->GetCell(0);
    FE3D FE_ID = fespace->GetFE3D(0, cell_0);
    TFE3D* FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
    RefTrans3D RefTrans = FE_Obj->GetRefTransID();

    // get the number of basis functions
  
}


void TParticles::CD_CC_Cuda()
{
    cout << "Inside CD_CC_Cuda" << endl;
    exit(0);
}

