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


__device__ void Obtain_velocity_at_a_point(int current_cell,
                                            int tid,
                                           double* d_m_particle_position_x,
                                           double* d_m_particle_position_y,
                                           double* d_m_particle_position_z,
                                           double* d_m_particle_velocity_x,
                                           double* d_m_particle_velocity_y,
                                           double* d_m_particle_velocity_z,
                                           double* d_m_cell_vertices_x,
                                           double* d_m_cell_vertices_y,
                                           double* d_m_cell_vertices_z,
                                            double* d_m_velocity_nodal_values_x,
                                            double* d_m_velocity_nodal_values_y,
                                            double* d_m_velocity_nodal_values_z,
                                           int* d_m_global_dof_indices,
                                           int* d_m_begin_indices)
{
        // -- Substep 1: Perform Setcell and get the reference values from original cells
        // -- Substep 2: Perform Interpolation of velocity values at the given particle position


        // -- Substep 1: Perform Setcell and get the reference values from original cells
        double x0, x1, x2, x3;
        double y0, y1, y2, y3;
        double z0, z1, z2, z3;

        // Fill the values of the vertices of the cell
        x0 = d_m_cell_vertices_x[4 * current_cell];
        x1 = d_m_cell_vertices_x[4 * current_cell + 1];
        x2 = d_m_cell_vertices_x[4 * current_cell + 2];
        x3 = d_m_cell_vertices_x[4 * current_cell + 3];

        y0 = d_m_cell_vertices_y[4 * current_cell];
        y1 = d_m_cell_vertices_y[4 * current_cell + 1];
        y2 = d_m_cell_vertices_y[4 * current_cell + 2];
        y3 = d_m_cell_vertices_y[4 * current_cell + 3];

        z0 = d_m_cell_vertices_z[4 * current_cell];
        z1 = d_m_cell_vertices_z[4 * current_cell + 1];
        z2 = d_m_cell_vertices_z[4 * current_cell + 2];
        z3 = d_m_cell_vertices_z[4 * current_cell + 3];


        double xc0=x0;
        double xc1=x1-x0;
        double xc2=x2-x0;
        double xc3=x3-x0;

        double yc0=y0;
        double yc1=y1-y0;
        double yc2=y2-y0;
        double yc3=y3-y0;

        double zc0=z0;
        double zc1=z1-z0;
        double zc2=z2-z0;
        double zc3=z3-z0;

        double detjk= xc1*yc2*zc3 + xc2*yc3*zc1 + xc3*yc1*zc2
                -xc3*yc2*zc1 - xc2*yc1*zc3 - xc1*yc3*zc2;
        
        // Sub block 2 : Get ref values from original cells
        double X = d_m_particle_position_x[tid];
        double Y = d_m_particle_position_y[tid];
        double Z = d_m_particle_position_z[tid];

        double xt=(X - xc0)/detjk;
        double yt=(Y - yc0)/detjk;
        double zt=(Z - zc0)/detjk;

        double xi  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
        double eta = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
        double zeta = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;

        // -- Get basis function values for the given cell
        double values[10];   // Hardcoded - THIVIN - the size of basis function is hardcoded for Tetraheadral p2
        
        
        double t1 = xi*xi;
        double t2 = xi*eta;
        double t3 = xi*zeta;
        double t4 = eta*eta;
        double t5 = eta*zeta;
        double t6 = zeta*zeta;

        values[0] = 1.0-3.0*xi-3.0*eta-3.0*zeta+2.0*t1+4.0*t2+4.0*t3
                    +2.0*t4+4.0*t5+2.0*t6;
        values[1] = 4.0*xi-4.0*t1-4.0*t2-4.0*t3;
        values[2] = -xi+2.0*t1;
        values[3] = 4.0*eta-4.0*t2-4.0*t4-4.0*t5;
        values[4] = 4.0*t2;
        values[5] = -eta+2.0*t4;
        values[6] = 4.0*zeta-4.0*t3-4.0*t5-4.0*t6;
        values[7] = 4.0*t3;
        values[8] = 4.0*t5;
        values[9] = -zeta+2.0*t6;


        // Lets get the start and end indices of the global dof indices
        int start_index = d_m_begin_indices[current_cell];
        int end_index = d_m_begin_indices[current_cell+1];
        
        // Now lets get the global dof indices
        int global_dof_indices[10];
        double nodal_velocity_x[10];
        double nodal_velocity_y[10];
        double nodal_velocity_z[10];

        double interpolated_velocity_x = 0.0;
        double interpolated_velocity_y = 0.0;
        double interpolated_velocity_z = 0.0;


        for(int i=0; i<10; i++)
        {
            global_dof_indices[i] = d_m_global_dof_indices[start_index+i];
            nodal_velocity_x[i] = d_m_velocity_nodal_values_x[global_dof_indices[i]];
            nodal_velocity_y[i] = d_m_velocity_nodal_values_y[global_dof_indices[i]];
            nodal_velocity_z[i] = d_m_velocity_nodal_values_z[global_dof_indices[i]];

            interpolated_velocity_x += values[i]*nodal_velocity_x[i];
            interpolated_velocity_y += values[i]*nodal_velocity_y[i];
            interpolated_velocity_z += values[i]*nodal_velocity_z[i];

            
        }

        

        // Assign the interpolated velocity to the particle
        d_m_particle_velocity_x[tid] = interpolated_velocity_x;
        d_m_particle_velocity_y[tid] = interpolated_velocity_y;
        d_m_particle_velocity_z[tid] = interpolated_velocity_z;

        // print the values
        printf("Particle GPU kernel 2-%d : %f %f %f\n", tid, interpolated_velocity_x, interpolated_velocity_y, interpolated_velocity_z);


}


__global__ void Interpolate_Velocity_CUDA(double* d_m_particle_position_x,
                                           double* d_m_particle_position_y,
                                           double* d_m_particle_position_z,
                                           double* d_m_particle_velocity_x,
                                           double* d_m_particle_velocity_y,
                                           double* d_m_particle_velocity_z,
                                           double* d_m_cell_vertices_x,
                                           double* d_m_cell_vertices_y,
                                           double* d_m_cell_vertices_z,
                                           double* d_m_velocity_nodal_values_x,
                                             double* d_m_velocity_nodal_values_y,
                                                double* d_m_velocity_nodal_values_z,
                                             int* d_m_current_cell,
                                                int* d_m_previous_cell,
                                           int* d_m_global_dof_indices,
                                           int* d_m_begin_indices,
                                           int n_cells,
                                           int n_dOF,
                                           int n_particles_released,
                                           double time_step)
{
    // Kernel code goes here
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if(tid < n_particles_released)
    {

        // -- BLOCK 1 : INterpolate velocity values at given particle position 
        // -- Substep 1: Identify the cell in which the particle is present


        // Get the cell_id_info
        double cell_no = d_m_current_cell[tid];
        double x0 = d_m_particle_position_x[tid];
        double y0 = d_m_particle_position_y[tid];
        double z0 = d_m_particle_position_z[tid];

        if (tid == 0)
        {
            // print first 10 velocity values
            for(int i=0; i<10; i++)
            {
                printf("Particle GPU kernel-1 : %d : %f %f %f\n", i, d_m_velocity_nodal_values_x[i], d_m_velocity_nodal_values_y[i], d_m_velocity_nodal_values_z[i]);
            }
        }

        // Call the funciton to obtain the velocity at the given point
        Obtain_velocity_at_a_point(cell_no,
                                    tid,
                                           d_m_particle_position_x,
                                           d_m_particle_position_y,
                                           d_m_particle_position_z,
                                           d_m_particle_velocity_x,
                                           d_m_particle_velocity_y,
                                           d_m_particle_velocity_z,
                                           d_m_cell_vertices_x,
                                           d_m_cell_vertices_y,
                                           d_m_cell_vertices_z,
                                           d_m_velocity_nodal_values_x,
                                             d_m_velocity_nodal_values_y,
                                                d_m_velocity_nodal_values_z,
                                           d_m_global_dof_indices,
                                           d_m_begin_indices);
    }

}

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
    h_m_cell_vertices_x = new double[4 * N_Cells];
    h_m_cell_vertices_y = new double[4 * N_Cells];
    h_m_cell_vertices_z = new double[4 * N_Cells];

    // Allocate memory for global indices of the cells
    h_m_global_dof_indices = new int[size_of_global_numbers];

    // Allocate memory for the begin index of the cells
    h_m_begin_indices  = new int[N_Cells + 1];

    // allocate memory for velocity values
    h_m_velocityX = new double[N_DOF];
    h_m_velocityY = new double[N_DOF];
    h_m_velocityZ = new double[N_DOF];


    // -- Now lets fill these values, so that these can be copied to the GPU

    // fill the cell vertices
    for(int i = 0; i < N_Cells; i++){
        // get the cell
        TBaseCell *cell = coll->GetCell(i);
        
        // fill the vertices
        for(int j = 0; j < 4; j++){
            double x0,y0,z0;
            cell->GetVertex(j)->GetCoords(x0, y0, z0);
            h_m_cell_vertices_x[4 * i + j] = x0;
            h_m_cell_vertices_y[4 * i + j] = y0;
            h_m_cell_vertices_z[4 * i + j] = z0;
        }
    }

    // Lets fill the global indices and begin index
    for(int i = 0; i < size_of_global_numbers; i++){
        h_m_global_dof_indices[i] = GlobalNumbers[i];
    }

    for(int i = 0; i < N_Cells + 1; i++){
        h_m_begin_indices[i] = BeginIndex[i];
    }

    // ----- ALLOCATE MEMORY IN GPU FOR ALL THE DEVICE VARIABLES -----

    // Allocate memory for the cell vertices
    checkCudaErrors(cudaMalloc((void**)&d_m_cell_vertices_x, 4 * N_Cells * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_cell_vertices_y, 4 * N_Cells * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_cell_vertices_z, 4 * N_Cells * sizeof(double)));

    // Allocate memory for global indices of the cells
    checkCudaErrors(cudaMalloc((void**)&d_m_global_dof_indices, size_of_global_numbers * sizeof(int)));

    // Allocate memory for the begin index of the cells
    checkCudaErrors(cudaMalloc((void**)&d_m_begin_indices, (N_Cells + 1) * sizeof(int)));

    // allocate memory for velocity values
    checkCudaErrors(cudaMalloc((void**)&d_m_velocity_nodal_values_x, N_DOF * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_velocity_nodal_values_y, N_DOF * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_velocity_nodal_values_z, N_DOF * sizeof(double)));

    // Allocate memory for current cell and previous cell
    checkCudaErrors(cudaMalloc((void**)&d_m_current_cell, N_Particles * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_m_previous_cell, N_Particles * sizeof(int)));

    // Allocate memory for the particle position
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_position_x, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_position_y, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_position_z, N_Particles * sizeof(double)));

    // Allocate memory for the particle previous position
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_previous_position_x, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_previous_position_y, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_previous_position_z, N_Particles * sizeof(double)));

    // Allocate memory for the particle velocity
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_velocity_x, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_velocity_y, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_velocity_z, N_Particles * sizeof(double)));


    // -- Copy the data from host to device -- //

    // copy the cell vertices
    checkCudaErrors(cudaMemcpy(d_m_cell_vertices_x, h_m_cell_vertices_x, 4 * N_Cells * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_cell_vertices_y, h_m_cell_vertices_y, 4 * N_Cells * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_cell_vertices_z, h_m_cell_vertices_z, 4 * N_Cells * sizeof(double), cudaMemcpyHostToDevice));

    // copy the global indices and begin index
    checkCudaErrors(cudaMemcpy(d_m_global_dof_indices, h_m_global_dof_indices, size_of_global_numbers * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_begin_indices, h_m_begin_indices, (N_Cells + 1) * sizeof(int), cudaMemcpyHostToDevice));

    // Copy the Initial particle position ( this is copied directly from the class variable, no parallel host array is created)
    checkCudaErrors(cudaMemcpy(d_m_particle_position_x, position_X.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_position_y, position_Y.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_position_z, position_Z.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));


    // Copy the previous particle position
    checkCudaErrors(cudaMemcpy(d_m_particle_previous_position_x, position_X_old.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_previous_position_y, position_Y_old.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_previous_position_z, position_Z_old.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));


    // Copy the particle velocity ( This will copy only zero values ), actual values will be copied in the main loop at each time step
    checkCudaErrors(cudaMemcpy(d_m_particle_velocity_x, velocityX.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_velocity_y, velocityY.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_velocity_z, velocityZ.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));

    // Copy the current cell and previous cell
    checkCudaErrors(cudaMemcpy(d_m_current_cell, currentCell.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_previous_cell, previousCell.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));


    cout << "[INFORMATION] Memory allocation and data transfer to GPU is done" << endl;
  
}


// Setup function to transfer velocity data at every time step
void TParticles::SetupVelocityValues(double *velocity_x_data, double *velocity_y_data, double *velocity_z_data, int N_particles_released, int N_DOF)
{
    // Copy the velocity values
    checkCudaErrors(cudaMemcpy(d_m_velocity_nodal_values_x, velocity_x_data, N_DOF * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_velocity_nodal_values_y, velocity_y_data, N_DOF * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_velocity_nodal_values_z, velocity_z_data, N_DOF * sizeof(double), cudaMemcpyHostToDevice));
}


void TParticles::CD_CC_Cuda()
{
    cout << "Inside CD_CC_Cuda" << endl;
    exit(0);
}


// Host wrapper for performing the velocity interpolation at every time step
void TParticles::InterpolateVelocityHostWrapper(double timeStep,int N_Particles_released,int N_DOF,int N_Cells)
{
    int MAX_THREAD_PER_BLOCK = 64;
    int N_threads; 

    if(N_Particles_released >= MAX_THREAD_PER_BLOCK) 
        N_threads = MAX_THREAD_PER_BLOCK;
        
    else
        N_threads = N_Particles_released;
    
    int C_NUM_BLOCKS = std::ceil(double(N_Particles_released)/MAX_THREAD_PER_BLOCK);


   dim3 dimGrid(C_NUM_BLOCKS);
   dim3 dimBlock(N_threads);

    Interpolate_Velocity_CUDA<<<dimGrid,dimBlock>>>(d_m_particle_position_x,
                                                     d_m_particle_position_y,
                                                     d_m_particle_position_z,
                                                     d_m_particle_velocity_x,
                                                     d_m_particle_velocity_y,
                                                     d_m_particle_velocity_z,
                                                     d_m_cell_vertices_x,
                                                     d_m_cell_vertices_y,
                                                     d_m_cell_vertices_z,
                                                     d_m_velocity_nodal_values_x,
                                                        d_m_velocity_nodal_values_y,
                                                        d_m_velocity_nodal_values_z,
                                                     d_m_current_cell,
                                                        d_m_previous_cell,
                                                     d_m_global_dof_indices,
                                                     d_m_begin_indices,
                                                     N_Cells,
                                                     N_DOF,
                                                     N_Particles_released,
                                                     timeStep);

    cudaDeviceSynchronize();

    // print the frst 10 values of the velocity for all 3 directions
    
    for(int i = 0; i < 10; i++)
    {
        cout <<"FINALCHECK : " <<velocityX[i] << " " << velocityY[i] << " " << velocityZ[i] << endl;
    }

    // Lets transfer the velocity values back to the host
    checkCudaErrors(cudaMemcpy(velocityX.data(), d_m_particle_velocity_x, N_Particles_released * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(velocityY.data(), d_m_particle_velocity_y, N_Particles_released * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(velocityZ.data(), d_m_particle_velocity_z, N_Particles_released * sizeof(double), cudaMemcpyDeviceToHost));

}




