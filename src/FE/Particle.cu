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




struct Vector3D_Cuda
{
    double x;
    double y;
    double z;
};





__device__ struct Vector3D_Cuda Obtain_velocity_at_a_point(int current_cell,
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

        // Create a new struct to store the interpolated velocity values
        struct Vector3D_Cuda interpolated_velocities;    

        // Assign the interpolated velocity to the particle
        interpolated_velocities.x = interpolated_velocity_x;
        interpolated_velocities.y = interpolated_velocity_y;
        interpolated_velocities.z = interpolated_velocity_z;

        // print the values
        // printf("Particle GPU kernel 2-%d : %f %f %f\n", tid, interpolated_velocity_x, interpolated_velocity_y, interpolated_velocity_z);
        return interpolated_velocities;
}

// GPU Function which computes the cd_cc value
__device__ double cd_cc_cuda(double fluid_density,
                          double particle_diameter,
                          double fluid_velocity,
                          double particle_velocity,
                          double fluid_dynamic_viscosity,
                          double lambda
                          )
{
    double Re_Particle = fluid_density * particle_diameter * abs(fluid_velocity - particle_velocity) / fluid_dynamic_viscosity;
    double cd = (24.0 / Re_Particle) * (1.0 + 0.15 * pow(Re_Particle, 0.687)); 
    double cc = 1.0 + ((2 * lambda) / particle_diameter) * (1.257 + 0.4 * exp(-1.0 * ((1.1 * particle_diameter) / (2 * lambda))));

    return cd/cc;
}   

// Function to set the sell and return the reference values
// this function is HARD CODED for a tetrahedral cell with p2 finite element
// HARDCODED - THIVIN
__device__ struct Vector3D_Cuda SetCell_And_Return_Reference_Value_CUDA(int current_cell,
                                                                    double* d_m_cell_vertices_x,
                                                                    double* d_m_cell_vertices_y,
                                                                    double* d_m_cell_vertices_z,
                                                                    double X,
                                                                    double Y,
                                                                    double Z)
{
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

        double xt=(X - xc0)/detjk;
        double yt=(Y - yc0)/detjk;
        double zt=(Z - zc0)/detjk;

        double xi  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
        double eta = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
        double zeta = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;

        // Create a new struct to store the xi eta zeta values
        struct Vector3D_Cuda xi_eta_zeta;

        // Assign the values
        xi_eta_zeta.x = xi;
        xi_eta_zeta.y = eta;
        xi_eta_zeta.z = zeta;

        return xi_eta_zeta;
}

// GPU function to calculate,if a point is within the given cell
// This routine is only for a tetrahedral cell with p2 Finite element
// HARDCODED - THIVIN
__device__ bool Is_Point_In_Cell_CUDA(int cellNo, 
                                    double* d_m_cell_vertices_x,
                                    double* d_m_cell_vertices_y,
                                    double* d_m_cell_vertices_z,
                                    double x,
                                    double y,
                                    double z)
{
    double xmin = 1e+8,  ymin = 1e+8, zmin = 1e+8;
    double xmax = -1e+8,  ymax = -1e+8, zmax = -1e+8;
    int i;
    bool ret = FALSE;

    double xi, eta, zeta;

    // Set cell Routine and get the reference values
    struct Vector3D_Cuda xi_eta_zeta = SetCell_And_Return_Reference_Value_CUDA(cellNo,
                                                                    d_m_cell_vertices_x,
                                                                    d_m_cell_vertices_y,
                                                                    d_m_cell_vertices_z,
                                                                    x,
                                                                    y,
                                                                    z);
    
    // Parse the values from the struct
    xi = xi_eta_zeta.x;
    eta = xi_eta_zeta.y;
    zeta = xi_eta_zeta.z;

    // Check if the reference values are within the tolerance range
    if(-1e-4 < xi && xi < 1.0001 &&
       -1e-4 < eta && eta < 1.0001 &&
       -1e-4 < zeta && zeta < 1.0001 &&
       xi + eta + zeta < 1.0001)
    {
      ret = TRUE;
    }

    return ret;
}


__global__ void Interpolate_Velocity_CUDA(  // cell Vertices
                                            double* d_m_cell_vertices_x,
                                            double* d_m_cell_vertices_y,
                                            double* d_m_cell_vertices_z,
                                            // Velocity Nodal Values
                                            double* d_m_velocity_nodal_values_x,
                                            double* d_m_velocity_nodal_values_y,
                                            double* d_m_velocity_nodal_values_z,
                                            // Particle Positions, Current & previous
                                            double* d_m_particle_position_x,
                                            double* d_m_particle_position_y,
                                            double* d_m_particle_position_z,
                                            double* d_m_particle_previous_position_x,
                                            double* d_m_particle_previous_position_y,
                                            double* d_m_particle_previous_position_z,
                                            // Particle Velocities, Current and Previous
                                            double* d_m_particle_velocity_x,
                                            double* d_m_particle_velocity_y,
                                            double* d_m_particle_velocity_z,
                                            double* d_m_particle_previous_velocity_x,
                                            double* d_m_particle_previous_velocity_y,
                                            double* d_m_particle_previous_velocity_z,
                                            // Particle & fluid Parameters
                                            double* d_m_particle_density,
                                            double* d_m_particle_diameter,
                                            double* d_m_fluid_density,
                                            double* d_m_dynamic_viscosity_fluid,
                                            double* d_m_lambda,
                                            double* d_m_gravity_x,
                                            double* d_m_gravity_y,
                                            double* d_m_gravity_z,
                                            //Computational Parameters
                                            int* d_m_current_cell,
                                            int* d_m_previous_cell,
                                            int* d_m_global_dof_indices,
                                            int* d_m_begin_indices,

                                            // FEM parameters for deposition
                                            int* d_m_is_boundary_cell,
                                            int* d_m_corner_id,
                                            int* d_m_is_boundary_dof_present,
                                            double* d_m_boundary_dof_x,
                                            double* d_m_boundary_dof_y,
                                            double* d_m_boundary_dof_z,
                                            int* d_m_joint_id,
                                            double* d_m_joint_normal_x,
                                            double* d_m_joint_normal_y,
                                            double* d_m_joint_normal_z,
                                            double* d_m_joint_coordinate_x,
                                            double* d_m_joint_coordinate_y,
                                            double* d_m_joint_coordinate_z,

                                            // STATISTICAL Parameters
                                            int* d_m_is_deposited_particle,
                                            int* d_m_is_escaped_particle,
                                            int* d_m_is_error_particle,
                                            int* d_m_is_stagnant_particle,
                                            int* d_m_is_ghost_particle,

                                            int n_cells,
                                            int n_dOF,
                                            int n_particles_released,
                                            double time_step
                                            )
{
    // Kernel code goes here
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if(tid < n_particles_released)
    {
        // -- BLOCK 1 : INterpolate velocity values at given particle position 
        // -- Substep 1: Identify the cell in which the particle is present

        // Get the cell_id_info
        int cell_no = d_m_current_cell[tid];
        double x0 = d_m_particle_position_x[tid];
        double y0 = d_m_particle_position_y[tid];
        double z0 = d_m_particle_position_z[tid];

        // if (tid == 0)
        // {
        //     printf("GPU : tid:  %d , cell_no: %d, x0:  %f, y0:  %f, z0: %f \n", tid, cell_no, x0, y0, z0);
        // }

        // Call the funciton to obtain the velocity at the given point
        struct Vector3D_Cuda interpolated_velocities = Obtain_velocity_at_a_point(cell_no,
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

        // USe the interpolated velocity functions to call the Calculate_Updated_Position_CUDA function 
        // This function updates the position array automatically
        // a Double can only be declated as a double * in cuda, due to malloc, so 
        // we will convert it to double value here.
        double fluid_density = *d_m_fluid_density;
        double dynamic_viscosity_fluid = *d_m_dynamic_viscosity_fluid;
        double lambda = *d_m_lambda;
        double gravity_x = *d_m_gravity_x;
        double gravity_y = *d_m_gravity_y;
        double gravity_z = *d_m_gravity_z;

        //  if (tid == 0)
        // {
        //     printf("GPU : %f %f %f\n", interpolated_velocities.x, interpolated_velocities.y, interpolated_velocities.z);
        // }

        // Calculate cd_cc
        double cd_cc_x = 0.0;
        double cd_cc_y = 0.0;
        double cd_cc_z = 0.0;

        // Assign te interpolated velocity values 
        double fluid_velocity_x = interpolated_velocities.x;
        double fluid_velocity_y = interpolated_velocities.y;
        double fluid_velocity_z = interpolated_velocities.z;


        cd_cc_x = cd_cc_cuda(fluid_density,
                            d_m_particle_diameter[tid],
                            fluid_velocity_x,
                            d_m_particle_velocity_x[tid],
                            dynamic_viscosity_fluid,
                            lambda
                            );
        
        cd_cc_y = cd_cc_cuda(fluid_density,
                                d_m_particle_diameter[tid],
                                fluid_velocity_y,
                                d_m_particle_velocity_y[tid],
                                dynamic_viscosity_fluid,
                                lambda
                                );
        
        cd_cc_z = cd_cc_cuda(fluid_density,
                                d_m_particle_diameter[tid],
                                fluid_velocity_z,
                                d_m_particle_velocity_z[tid],
                                dynamic_viscosity_fluid,
                                lambda
                                );  
        
        // FInd the minimum cd_cc
        double cd_cc_min = cd_cc_x;
        if(cd_cc_y < cd_cc_min)
            cd_cc_min = cd_cc_y;
        if(cd_cc_z < cd_cc_min)
            cd_cc_min = cd_cc_z;
        
        // Set the cd_cc_x, cd_cc_y, cd_cc_z to cd_cc_min
        cd_cc_x = cd_cc_min;
        cd_cc_y = cd_cc_min;
        cd_cc_z = cd_cc_min;

        // Calculate the RHS 
        double rhs_x = 0.0;
        double rhs_y = 0.0;
        double rhs_z = 0.0;

        double inertial_constant = (3. / 4.) * (fluid_density / d_m_particle_density[tid]) * (1 / d_m_particle_diameter[tid]);

        rhs_x = inertial_constant * cd_cc_x * abs(fluid_velocity_x - d_m_particle_velocity_x[tid]) * (fluid_velocity_x - d_m_particle_velocity_x[tid]);
        rhs_x += gravity_x * (fluid_density - d_m_particle_density[tid]) / d_m_particle_density[tid];

        rhs_y = inertial_constant * cd_cc_y * abs(fluid_velocity_y - d_m_particle_velocity_y[tid]) * (fluid_velocity_y - d_m_particle_velocity_y[tid]);
        rhs_y += gravity_y * (fluid_density - d_m_particle_density[tid]) / d_m_particle_density[tid];

        rhs_z = inertial_constant * cd_cc_z * abs(fluid_velocity_z - d_m_particle_velocity_z[tid]) * (fluid_velocity_z - d_m_particle_velocity_z[tid]);
        rhs_z += gravity_z * (fluid_density - d_m_particle_density[tid]) / d_m_particle_density[tid];

        //  if (tid == 0)
        // {
        //     printf("GPU : %f %f %f\n", rhs_x, rhs_y, rhs_z);
        // }

        // Compute the updated paticle velocity using forward euler
        d_m_particle_velocity_x[tid] = rhs_x * time_step + d_m_particle_velocity_x[tid];
        d_m_particle_velocity_y[tid] = rhs_y * time_step + d_m_particle_velocity_y[tid];
        d_m_particle_velocity_z[tid] = rhs_z * time_step + d_m_particle_velocity_z[tid];

        // if (tid == 0)
        // {
        //     printf("GPU : %f %f %f\n", d_m_particle_velocity_x[tid], d_m_particle_velocity_y[tid], d_m_particle_velocity_z[tid]);
        // }

        // Update the particle position using RK-2
        d_m_particle_position_x[tid] = time_step * 0.5 * (d_m_particle_velocity_x[tid] + d_m_particle_previous_velocity_x[tid]) + d_m_particle_position_x[tid];
        d_m_particle_position_y[tid] = time_step * 0.5 * (d_m_particle_velocity_y[tid] + d_m_particle_previous_velocity_y[tid]) + d_m_particle_position_y[tid];
        d_m_particle_position_z[tid] = time_step * 0.5 * (d_m_particle_velocity_z[tid] + d_m_particle_previous_velocity_z[tid]) + d_m_particle_position_z[tid];
        

        // Transfer the updated particle velocity to the previous particle velocity
        d_m_particle_previous_velocity_x[tid] = d_m_particle_velocity_x[tid];
        d_m_particle_previous_velocity_y[tid] = d_m_particle_velocity_y[tid];
        d_m_particle_previous_velocity_z[tid] = d_m_particle_velocity_z[tid];

        // if (tid == 0)
        // {
        //     printf("GPU : %f %f %f\n", d_m_particle_previous_velocity_x[tid], d_m_particle_previous_velocity_y[tid], d_m_particle_previous_velocity_z[tid]);
        // }

        // Transfer current particle position to previous particle position
        d_m_particle_previous_position_x[tid] = d_m_particle_position_x[tid];
        d_m_particle_previous_position_y[tid] = d_m_particle_position_y[tid];
        d_m_particle_previous_position_z[tid] = d_m_particle_position_z[tid];


        // Check the current position of the cells within the domain 
        bool inside_domain = false;
        for (int cell_id = 0; cell_id < n_cells; cell_id++)
        {
            bool insideCell = Is_Point_In_Cell_CUDA(cell_id,
                                                    d_m_cell_vertices_x,
                                                    d_m_cell_vertices_y,
                                                    d_m_cell_vertices_z,
                                                    d_m_particle_position_x[tid],
                                                    d_m_particle_position_y[tid],
                                                    d_m_particle_position_z[tid]);
            if (insideCell)
            {
                inside_domain = true;
                // copy the current cell to previous cell
                d_m_previous_cell[tid] = d_m_current_cell[tid];
                d_m_current_cell[tid] = cell_id;
                
                break;
            }
        }

        // If not inside the domain, then the particle is either escaped or deposited.
        if(!inside_domain)
        {
           // Find the last cell that the particle was present in 
           int cell_no = d_m_current_cell[tid];

            // Check if the cell is a boundary cell (cell with a boundary face)
            int is_boundary_cell = d_m_is_boundary_cell[cell_no];
            int index_boundary_cell = 0;
            int index_bounday_dof = 0;
            
            
            if(is_boundary_cell < 0) // non positive value, preferably -99999, which is used as place holder
                is_boundary_cell = 0;
            else
                index_boundary_cell = is_boundary_cell; // assign the index

            // Check if boundary DOF is present
            int is_boundary_dof_present = d_m_is_boundary_dof_present[tid];

            if(is_boundary_dof_present < 0) // non positive value, preferably -99999, which is used as place holder
                is_boundary_dof_present = 0;
            else
                index_bounday_dof = is_boundary_dof_present; // assign the index

            int corner_id;
            int joint_id;
            
            // If last cell not a boundary cell, check for the boundary DOF's within the cell. 
            if(!is_boundary_cell)
            {
                if(is_boundary_dof_present)  // Atleast a boundary DOF is present
                {
                    // obtain the bound dof coordinates using the index_bounday_dof index from d_m_boundary_dof_x, d_m_boundary_dof_y, d_m_boundary_dof_z
                    double boundary_dof_x = d_m_boundary_dof_x[index_bounday_dof];
                    double boundary_dof_y = d_m_boundary_dof_y[index_bounday_dof];
                    double boundary_dof_z = d_m_boundary_dof_z[index_bounday_dof];

                    // Assign it to the particle position
                    d_m_particle_position_x[tid] = boundary_dof_x;
                    d_m_particle_position_y[tid] = boundary_dof_y;
                    d_m_particle_position_z[tid] = boundary_dof_z;

                    // Mark the particle as deposited
                    d_m_is_deposited_particle[tid] = 1;
                    return;
                }

                // FUTURE TODO : May be add a check to see if the neighbouring cells have boundary DOF's or boundary faces
                else  // Does not have a Boundary face nor a boundary DOF
                {
                    // Mark the particle as escaped
                    d_m_is_deposited_particle[tid] = 1;
                    return;
                }
            }

            // if the corner id = 21, then it means, the cell has boundary faces on wall and outflow boundary
            // So the particle is at the the last cell near the exit
            // since the particle has reached here, we will consider the particle escaped via the outflow boundary
            if (is_boundary_cell && corner_id == 21 ) 
            {
                // Mark the particle as deposited , for book keeping purposes
                d_m_is_deposited_particle[tid] = 1;

                // Mark the particle as escaped
                d_m_is_escaped_particle[tid] = 1;

                return;
            }

             // Check if the particle is from a corner shared by two bdids 2 and 0
             // this is an inlet between wall and inlet, so we will consider the particle escaped via the wall
            if (is_boundary_cell && corner_id == 20)
            {
                joint_id  = d_m_joint_id[index_boundary_cell];  
            }

            // If a paricle escaped from a cell, in which there is only inlet surface then its an error particle
            // This will only happen if there is a backflow
            if(corner_id == 0)
            {
                // Mark the particle as deposited , for book keeping purposes
                d_m_is_deposited_particle[tid] = 1;
                
                // Mark the particle as error 
                d_m_is_error_particle[tid] = 1;
                
                return;
            }
            else // Corner ID is 1 or 2
            {   
                // get the joint id
                joint_id  = d_m_joint_id[index_boundary_cell];
            }

            // With the obtained joint id and get the corresponding joint normal and joint coordinate
            double joint_normal_x = d_m_joint_normal_x[index_boundary_cell];
            double joint_normal_y = d_m_joint_normal_y[index_boundary_cell];
            double joint_normal_z = d_m_joint_normal_z[index_boundary_cell];

            double joint_coordinate_x = d_m_joint_coordinate_x[index_boundary_cell];
            double joint_coordinate_y = d_m_joint_coordinate_y[index_boundary_cell];
            double joint_coordinate_z = d_m_joint_coordinate_z[index_boundary_cell];

            // Get the structs
            struct Vector3D_Cuda first_position;
            struct Vector3D_Cuda second_position;
            struct Vector3D_Cuda line_vector;
            struct Vector3D_Cuda point_on_surface;
            struct Vector3D_Cuda normal_surface;
            struct Vector3D_Cuda temp1;
            struct Vector3D_Cuda temp2;

            first_position.x = d_m_particle_previous_position_x[tid];
            first_position.y = d_m_particle_previous_position_y[tid];
            first_position.z = d_m_particle_previous_position_z[tid];

            second_position.x = d_m_particle_position_x[tid];
            second_position.y = d_m_particle_position_y[tid];
            second_position.z = d_m_particle_position_z[tid];

            normal_surface.x = joint_normal_x;
            normal_surface.y = joint_normal_y;
            normal_surface.z = joint_normal_z;

            point_on_surface.x = joint_coordinate_x;
            point_on_surface.y = joint_coordinate_y;
            point_on_surface.z = joint_coordinate_z;

            //  u = p1 - p0
            line_vector.x = first_position.x - second_position.x;
            line_vector.y = first_position.y - second_position.y;
            line_vector.z = first_position.z - second_position.z;

            // Dot
            double dot = normal_surface.x * first_position.x + normal_surface.y * first_position.y + normal_surface.z * first_position.z;

            if(fabs(dot - 0.0) > 1e-3)
            {
                // w = p0 - pC0
                temp1.x = first_position.x - point_on_surface.x;
                temp1.y = first_position.y - point_on_surface.y;
                temp1.z = first_position.z - point_on_surface.z;

                double fac = -1.0 * (normal_surface.x * temp1.x + normal_surface.y * temp1.y + normal_surface.z * temp1.z);
                fac /= fac;

                // u = u*fac
                line_vector.x = line_vector.x * fac;
                line_vector.y = line_vector.y * fac;
                line_vector.z = line_vector.z * fac;

                temp2.x = line_vector.x + first_position.x;
                temp2.y = line_vector.y + first_position.y;
                temp2.z = line_vector.z + first_position.z;

                // Mark the Particle as deposited
                d_m_is_deposited_particle[tid] = 1;

                d_m_particle_position_x[tid] = temp2.x;
                d_m_particle_position_y[tid] = temp2.y;
                d_m_particle_position_z[tid] = temp2.z;

                // if bdId/corner ID is 1, then mark the particle as escaped
                if(corner_id == 1)
                {
                    // Mark the Particle as deposited
                    d_m_is_deposited_particle[tid] = 1;

                    // Mark the particle as escaped
                    d_m_is_escaped_particle[tid] = 1;

                    return;
                }
            }
            else
            {
                if(corner_id ==1)
                {
                    d_m_is_deposited_particle[tid] = 1; // Mark the Particle as deposited
                    d_m_is_escaped_particle[tid] = 1;   // Mark the particle as escaped
                }

                // mark them as deposited 

                d_m_is_deposited_particle[tid] = 1;

                // Ghost particle, make the vertex position as deposition position
                d_m_particle_position_x[tid] = point_on_surface.x;
                d_m_particle_position_y[tid] = point_on_surface.y;
                d_m_particle_position_z[tid] = point_on_surface.z;

                d_m_is_ghost_particle[tid] = 1;
                
            }

        }

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

    // Refer to Documentation "DepositionLogic.md" for further details. 
    // In a Nutshell, there will be a `is_***` array, which will be of size N_Cells
    // if the cell does not match the category, then the value will be given as `-99999`
    // if the cell do match the category, then the value will be given as a positive number, which will be the index at which the corresponding values will be present

    // For Eg : Lets say , there are 5 cells , and cells 0,2,3 has boundary_faces, then the is_boundary_faces array will be as follows
    // is_boundary_faces = [0, -99999,1, 2, -99999], and the joint ids of the cells 0,2,3 wil be saved in the array
    // joint_ids = [4, 0, 1] , where 4 is the joint id of cell 0, 0 is the joint id of cell 2 and 1 is the joint id of cell 3

    // Lets fill the is_boundary_faces array
    h_m_is_boundary_cell = new int[N_Cells];

    // get size of the m_mapBoundaryFaceIds map
    int size_of_boundary_face = m_mapBoundaryFaceIds.size();
    
    // create joint_ids array and corner_ids array of size equal to the size of the m_mapBoundaryFaceIds map
    h_m_joint_id = new int[size_of_boundary_face];
    h_m_corner_id = new int[size_of_boundary_face];

    // Allocate memory for joint normals and joint coordinates
    h_m_joint_normal_x = new double[size_of_boundary_face];
    h_m_joint_normal_y = new double[size_of_boundary_face];
    h_m_joint_normal_z = new double[size_of_boundary_face];

    h_m_joint_coordinate_x = new double[size_of_boundary_face];
    h_m_joint_coordinate_y = new double[size_of_boundary_face];
    h_m_joint_coordinate_z = new double[size_of_boundary_face];


    // Initialise the h_m_is_boundary_faces, h_m_joint_ids and h_m_corner_ids array 
    // using maps, m_mapBoundaryFaceIds, m_mapJointIds and m_mapCornerIds
    // if the key is present in the map, then the value will be the index at which the corresponding values will be present
    // if the key is not present in the map, then the value will be -99999

    int index_for_boundary_face = 0;
    // loop through the cells
    for (int cell_no = 0 ; cell_no < N_Cells; cell_no++)
    {
        // check if the cell no is present in the map
        if (m_mapBoundaryFaceIds.find(cell_no) != m_mapBoundaryFaceIds.end())
        {
            // if the cell no is present in the map, then the value will be the index at which the corresponding values will be present
            h_m_is_boundary_cell[cell_no] = index_for_boundary_face;

            // get the joint id and corner id from the map
            int joint_id = m_jointidOfBoundCells[cell_no];
            int corner_id = m_cornerTypeOfBoundCells[cell_no];

            // save the joint id and corner id in the corresponding arrays
            h_m_joint_id[index_for_boundary_face] = joint_id;
            h_m_corner_id[index_for_boundary_face] = corner_id;

            // get the joint normal and joint coordinate from the map
            double joint_normal_x ;
            double joint_normal_y ;
            double joint_normal_z ;

            double joint_coordinate_x;
            double joint_coordinate_y;
            double joint_coordinate_z;


            // Variables required for normal Computation 
            int MaxLen;
            int N_Joints;
            const int *TmpLen;
            const int *TmpFV;


            TBaseCell *cell = fespace->GetCollection()->GetCell(cell_no);
            cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
            TJoint *Joint = cell->GetJoint(joint_id);
            double x1, x2, x3, y1, y2, y3, z1, z2, z3;

            // Get the coordinates of the joint
            cell->GetVertex(TmpFV[joint_id * MaxLen + 0])->GetCoords(x1, y1, z1);
            cell->GetVertex(TmpFV[joint_id * MaxLen + 1])->GetCoords(x2, y2, z2);
            double t11 = x2 - x1;
            double t12 = y2 - y1;
            double t13 = z2 - z1;
            double len = sqrt(t11 * t11 + t12 * t12 + t13 * t13);
            t11 /= len;
            t12 /= len;
            t13 /= len;

            cell->GetVertex(TmpFV[joint_id * MaxLen + (TmpLen[joint_id] - 1)])->GetCoords(x2, y2, z2);
            double t21 = x2 - x1;
            double t22 = y2 - y1;
            double t23 = z2 - z1;
            len = sqrt(t21 * t21 + t22 * t22 + t23 * t23);
            t21 /= len;
            t22 /= len;
            t23 /= len;

            double N1 = t12 * t23 - t13 * t22;
            double N2 = t13 * t21 - t11 * t23;
            double N3 = t11 * t22 - t12 * t21;
            len = sqrt(N1 * N1 + N2 * N2 + N3 * N3);
            N1 /= len;
            N2 /= len;
            N3 /= len;

            // Assign the joint normal 
            joint_normal_x = N1;
            joint_normal_y = N2;
            joint_normal_z = N3;

            // Assign joint coordinate
            joint_coordinate_x = x1;
            joint_coordinate_y = y1;
            joint_coordinate_z = z1;

            // save the joint normal and joint coordinate in the corresponding arrays
            h_m_joint_normal_x[index_for_boundary_face] = joint_normal_x;
            h_m_joint_normal_y[index_for_boundary_face] = joint_normal_y;
            h_m_joint_normal_z[index_for_boundary_face] = joint_normal_z;

            h_m_joint_coordinate_x[index_for_boundary_face] = joint_coordinate_x;
            h_m_joint_coordinate_y[index_for_boundary_face] = joint_coordinate_y;
            h_m_joint_coordinate_z[index_for_boundary_face] = joint_coordinate_z;

            // increment the index_for_boundary_face
            index_for_boundary_face++;
        }
        else
        {
            // if the cell no is not present in the map, then the value will be -99999
            h_m_is_boundary_cell[cell_no] = -99999;
        }
    }


    // fill the is_boundary_dof array 
    h_m_is_boundary_dof_present = new int[N_Cells];

    // get the size of the m_mapBoundaryDofIds map
    int size_of_boundary_dof = m_BoundaryDOFsOnCell.size();

    // Allocate memory for the arrays
    h_m_boundary_dof_x = new double[size_of_boundary_dof];
    h_m_boundary_dof_y = new double[size_of_boundary_dof];
    h_m_boundary_dof_z = new double[size_of_boundary_dof];
    
    // Initialise the h_m_is_boundary_dof, h_m_joint_normal_x, h_m_joint_normal_y, h_m_joint_normal_z, h_m_joint_coordinate_x, h_m_joint_coordinate_y, h_m_joint_coordinate_z array
    // using maps, m_mapBoundaryDofIds, and compute the joint normals and joint coordinates
    
    // print size of boundary dof
    cout << "[INFORMATION] size_of_boundary_dof : " << size_of_boundary_dof << endl;

    // print num of cells
    cout << "[INFORMATION] N_Cells : " << N_Cells << endl;


    int index_for_boundary_dof = 0;
    
    for (int cell_no = 0 ; cell_no < N_Cells ; cell_no ++)
    {
        // check if the cell no is present in the map
        if (m_BoundaryDOFsOnCell.find(cell_no) != m_BoundaryDOFsOnCell.end())
        {
            // if the cell no is present in the map, then the value will be the index at which the corresponding values will be present
            h_m_is_boundary_dof_present[cell_no] = index_for_boundary_dof;
            
            std::vector<double> boundary_dof = m_BoundaryDOFsOnCell[cell_no];
            h_m_boundary_dof_x[index_for_boundary_dof] = boundary_dof[0]; // x dof
            h_m_boundary_dof_y[index_for_boundary_dof] = boundary_dof[1]; // y dof
            h_m_boundary_dof_z[index_for_boundary_dof] = boundary_dof[2]; // z dof

            // increment the index_for_boundary_dof
            index_for_boundary_dof++;
        }
        else
        {
            // if the cell no is not present in the map, then the value will be -99999
            h_m_is_boundary_dof_present[cell_no] = -99999;
        }
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

    // Allocate memory for the particle previous velocity
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_previous_velocity_x, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_previous_velocity_y, N_Particles * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_previous_velocity_z, N_Particles * sizeof(double)));

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


    // Copy the previous particle velocity
    checkCudaErrors(cudaMemcpy(d_m_particle_previous_velocity_x, velocityX_old.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_previous_velocity_y, velocityY_old.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_particle_previous_velocity_z, velocityZ_old.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));

    // Copy the current cell and previous cell
    checkCudaErrors(cudaMemcpy(d_m_current_cell, currentCell.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_m_previous_cell, previousCell.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));


    // ----- PARTICLE AND FLUID PARAMETERS ----- //
    
    // Allocate memory for particle density
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_density, N_Particles * sizeof(double)));

    // Allocate memory for particle diameter
    checkCudaErrors(cudaMalloc((void**)&d_m_particle_diameter, N_Particles * sizeof(double)));

    // Allocate memory for fluid density : Data type is scalar and double
    checkCudaErrors(cudaMalloc((void**)&d_m_fluid_density, sizeof(double)));

    // Allocate memory for fluid viscosity : Data type is scalar and double
    checkCudaErrors(cudaMalloc((void**)&d_m_dynamic_viscosity_fluid, sizeof(double)));

    // Allocate memory for fluid gravity in x direction : Data type is scalar and double
    checkCudaErrors(cudaMalloc((void**)&d_m_gravity_x, sizeof(double)));

    // Allocate memory for fluid gravity in y direction : Data type is scalar and double
    checkCudaErrors(cudaMalloc((void**)&d_m_gravity_y, sizeof(double)));

    // Allocate memory for fluid gravity in z direction : Data type is scalar and double
    checkCudaErrors(cudaMalloc((void**)&d_m_gravity_z, sizeof(double)));

    // Allocate lambda : Data type is scalar and double
    checkCudaErrors(cudaMalloc((void**)&d_m_lambda, sizeof(double)));

    // Send the fluid and particle parameters to GPU, 
    // d_m_ is the device variable with dtype as double* and m_ is the host variable with dtype as std::vector<double> 
    checkCudaErrors(cudaMemcpy(d_m_particle_density, m_particle_density.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));

    // Send the particle diameter
    checkCudaErrors(cudaMemcpy(d_m_particle_diameter, m_particle_diameter.data(), N_Particles * sizeof(double), cudaMemcpyHostToDevice));

    // Send the fluid density
    checkCudaErrors(cudaMemcpy(d_m_fluid_density, &m_fluid_density, sizeof(double), cudaMemcpyHostToDevice));

    // Send the fluid viscosity
    checkCudaErrors(cudaMemcpy(d_m_dynamic_viscosity_fluid, &m_fluid_dynamic_viscosity, sizeof(double), cudaMemcpyHostToDevice));

    // Send the fluid gravity in x direction
    checkCudaErrors(cudaMemcpy(d_m_gravity_x, &m_gravity_x, sizeof(double), cudaMemcpyHostToDevice));

    // Send the fluid gravity in y direction
    checkCudaErrors(cudaMemcpy(d_m_gravity_y, &m_gravity_y, sizeof(double), cudaMemcpyHostToDevice));

    // Send the fluid gravity in z direction
    checkCudaErrors(cudaMemcpy(d_m_gravity_z, &m_gravity_z, sizeof(double), cudaMemcpyHostToDevice));

    // Send the lambda
    checkCudaErrors(cudaMemcpy(d_m_lambda, &m_lambda, sizeof(double), cudaMemcpyHostToDevice));

    // -- TIME PARAMETERS -- //

    // Allocate memory for the time step
    checkCudaErrors(cudaMalloc((void**)&d_m_time_step, sizeof(double)));

    // transfer the time step to the GPU
    checkCudaErrors(cudaMemcpy(d_m_time_step, &h_m_time_step, sizeof(double), cudaMemcpyHostToDevice));


    // ---- DEPOSITION PARAMETERS -- // 
    // Allocate memory for d_m_is_boundary_cell with size N_Cells and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_boundary_cell, N_Cells * sizeof(int)));

    // Allocate memory for d_m_corner_id with size equal to the size of the m_mapBoundaryFaceIds map and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_corner_id, size_of_boundary_face * sizeof(int)));

    // Allocate memory for d_m_joint_id with size equal to the size of the m_mapBoundaryFaceIds map and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_id, size_of_boundary_face * sizeof(int)));

    // Allocate memory for d_m_joint_normal_x with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_normal_x, size_of_boundary_face * sizeof(double)));

    // Allocate memory for d_m_joint_normal_y with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_normal_y, size_of_boundary_face * sizeof(double)));

    // Allocate memory for d_m_joint_normal_z with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_normal_z, size_of_boundary_face * sizeof(double)));

    // Allocate memory for d_m_joint_coordinate_x with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_coordinate_x, size_of_boundary_face * sizeof(double)));

    // Allocate memory for d_m_joint_coordinate_y with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_coordinate_y, size_of_boundary_face * sizeof(double)));

    // Allocate memory for d_m_joint_coordinate_z with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMalloc((void**)&d_m_joint_coordinate_z, size_of_boundary_face * sizeof(double)));

    // Allocate memory for d_m_is_boundary_dof_present with size N_cells and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_boundary_dof_present, N_Cells * sizeof(int)));

    // Allocate memory for d_m_boundary_dof_x with size equal to the size of the m_mapBoundaryDofIds map and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_boundary_dof_x, size_of_boundary_dof * sizeof(int)));

    // Allocate memory for d_m_boundary_dof_y with size equal to the size of the m_mapBoundaryDofIds map and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_boundary_dof_y, size_of_boundary_dof * sizeof(int)));

    // Allocate memory for d_m_boundary_dof_z with size equal to the size of the m_mapBoundaryDofIds map and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_boundary_dof_z, size_of_boundary_dof * sizeof(int)));
    

    // Copy the values from host to these arrays 
    // copy for d_m_is_boundary_cell from h_m_is_boundary_cell with size N_Cells and data type int
    checkCudaErrors(cudaMemcpy(d_m_is_boundary_cell, h_m_is_boundary_cell, N_Cells * sizeof(int), cudaMemcpyHostToDevice));

    // copy for d_m_corner_id from h_m_corner_id with size equal to the size of the m_mapBoundaryFaceIds map and data type int
    checkCudaErrors(cudaMemcpy(d_m_corner_id, h_m_corner_id, size_of_boundary_face * sizeof(int), cudaMemcpyHostToDevice));

    // copy for d_m_joint_id from h_m_joint_id with size equal to the size of the m_mapBoundaryFaceIds map and data type int
    checkCudaErrors(cudaMemcpy(d_m_joint_id, h_m_joint_id, size_of_boundary_face * sizeof(int), cudaMemcpyHostToDevice));

    // copy for d_m_joint_normal_x from h_m_joint_normal_x with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMemcpy(d_m_joint_normal_x, h_m_joint_normal_x, size_of_boundary_face * sizeof(double), cudaMemcpyHostToDevice));

    // copy for d_m_joint_normal_y from h_m_joint_normal_y with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMemcpy(d_m_joint_normal_y, h_m_joint_normal_y, size_of_boundary_face * sizeof(double), cudaMemcpyHostToDevice));

    // copy for d_m_joint_normal_z from h_m_joint_normal_z with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMemcpy(d_m_joint_normal_z, h_m_joint_normal_z, size_of_boundary_face * sizeof(double), cudaMemcpyHostToDevice));

    // copy for d_m_joint_coordinate_x from h_m_joint_coordinate_x with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMemcpy(d_m_joint_coordinate_x, h_m_joint_coordinate_x, size_of_boundary_face * sizeof(double), cudaMemcpyHostToDevice));

    // copy for d_m_joint_coordinate_y from h_m_joint_coordinate_y with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMemcpy(d_m_joint_coordinate_y, h_m_joint_coordinate_y, size_of_boundary_face * sizeof(double), cudaMemcpyHostToDevice));

    // copy for d_m_joint_coordinate_z from h_m_joint_coordinate_z with size equal to the size of the m_mapBoundaryFaceIds map and data type double
    checkCudaErrors(cudaMemcpy(d_m_joint_coordinate_z, h_m_joint_coordinate_z, size_of_boundary_face * sizeof(double), cudaMemcpyHostToDevice));

    // copy for d_m_is_boundary_dof_present from h_m_is_boundary_dof_present with size N_cells and data type int
    checkCudaErrors(cudaMemcpy(d_m_is_boundary_dof_present, h_m_is_boundary_dof_present, N_Cells * sizeof(int), cudaMemcpyHostToDevice));

    // copy for d_m_boundary_dof_x from h_m_boundary_dof_x with size equal to the size of the m_mapBoundaryDofIds map and data type int
    checkCudaErrors(cudaMemcpy(d_m_boundary_dof_x, h_m_boundary_dof_x, size_of_boundary_dof * sizeof(int), cudaMemcpyHostToDevice));

    // copy for d_m_boundary_dof_y from h_m_boundary_dof_y with size equal to the size of the m_mapBoundaryDofIds map and data type int
    checkCudaErrors(cudaMemcpy(d_m_boundary_dof_y, h_m_boundary_dof_y, size_of_boundary_dof * sizeof(int), cudaMemcpyHostToDevice));

    // copy for d_m_boundary_dof_z from h_m_boundary_dof_z with size equal to the size of the m_mapBoundaryDofIds map and data type int
    checkCudaErrors(cudaMemcpy(d_m_boundary_dof_z, h_m_boundary_dof_z, size_of_boundary_dof * sizeof(int), cudaMemcpyHostToDevice));

    // --  PARTICLE LEVEL STATISTICAL VARIABLES -- //

    // Allocate memory for the d_m_is_escaped_particle array, size N_Particles and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_escaped_particle, N_Particles * sizeof(int)));
    
    // Allocate memory for the d_m_is_deposited_particle array, size N_Particles and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_deposited_particle, N_Particles * sizeof(int)));

    // Allocate memory for the d_m_is_error_particle array, size N_Particles and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_error_particle, N_Particles * sizeof(int)));

    // Allocate memory for the d_m_is_stagnant_particle array, size N_Particles and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_stagnant_particle, N_Particles * sizeof(int)));

    // Allocate memory for the d_m_ghost_particle array, size N_Particles and data type int
    checkCudaErrors(cudaMalloc((void**)&d_m_is_ghost_particle, N_Particles * sizeof(int)));


    // Copy the values from the host to the device

    // copy the isParticleDeposited vector to d_m_is_deposited_particle with size N_Particles and data type int
    checkCudaErrors(cudaMemcpy(d_m_is_deposited_particle, isParticleDeposited.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));

    // copy the isEscapedParticle vector to d_m_is_escaped_particle with size N_Particles and data type int
    checkCudaErrors(cudaMemcpy(d_m_is_escaped_particle, isEscapedParticle.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));

    // copy the isStagnantParticle vector to d_m_is_stagnant_particle with size N_Particles and data type int
    checkCudaErrors(cudaMemcpy(d_m_is_stagnant_particle, isStagnantParticle.data(), N_Particles * sizeof(int), cudaMemcpyHostToDevice));

    // initialize the d_m_is_ghost_particle array with 0 on the device  // There is no Ghost particle array in the given CPU code
    checkCudaErrors(cudaMemset(d_m_is_ghost_particle, 0, N_Particles * sizeof(int)));

    

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
void TParticles::InterpolateVelocityHostWrapper(double time_step,int N_Particles_released,int N_DOF,int N_Cells)
{
    int MAX_THREAD_PER_BLOCK = 128;
    int N_threads; 

    if(N_Particles_released >= MAX_THREAD_PER_BLOCK) 
        N_threads = MAX_THREAD_PER_BLOCK;
        
    else
        N_threads = N_Particles_released;
    
    int C_NUM_BLOCKS = std::ceil(double(N_Particles_released)/MAX_THREAD_PER_BLOCK);


   dim3 dimGrid(C_NUM_BLOCKS);
   dim3 dimBlock(N_threads);

   // time the kernel
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventRecord(start, 0);

    Interpolate_Velocity_CUDA<<<dimGrid,dimBlock>>>(// cell Vertices
                                                    d_m_cell_vertices_x,
                                                    d_m_cell_vertices_y,
                                                    d_m_cell_vertices_z,
                                                    // Velocity Nodal Values
                                                    d_m_velocity_nodal_values_x,
                                                    d_m_velocity_nodal_values_y,
                                                    d_m_velocity_nodal_values_z,
                                                    // Particle Positions, Current & previous
                                                    d_m_particle_position_x,
                                                    d_m_particle_position_y,
                                                    d_m_particle_position_z,
                                                    d_m_particle_previous_position_x,
                                                    d_m_particle_previous_position_y,
                                                    d_m_particle_previous_position_z,
                                                    // Particle Velocities, Current and Previous
                                                    d_m_particle_velocity_x,
                                                    d_m_particle_velocity_y,
                                                    d_m_particle_velocity_z,
                                                    d_m_particle_previous_velocity_x,
                                                    d_m_particle_previous_velocity_y,
                                                    d_m_particle_previous_velocity_z,
                                                    // Particle & fluid Parameters
                                                    d_m_particle_density,
                                                    d_m_particle_diameter,
                                                    d_m_fluid_density,
                                                    d_m_dynamic_viscosity_fluid,
                                                    d_m_lambda,
                                                    d_m_gravity_x,
                                                    d_m_gravity_y,
                                                    d_m_gravity_z,
                                                    //Computational Parameters
                                                    d_m_current_cell,
                                                    d_m_previous_cell,
                                                    d_m_global_dof_indices,
                                                    d_m_begin_indices,

                                                    // FEM parameters for deposition
                                                    d_m_is_boundary_cell,
                                                    d_m_corner_id,
                                                    d_m_is_boundary_dof_present,
                                                    d_m_boundary_dof_x,
                                                    d_m_boundary_dof_y,
                                                    d_m_boundary_dof_z,
                                                    d_m_joint_id,
                                                    d_m_joint_normal_x,
                                                    d_m_joint_normal_y,
                                                    d_m_joint_normal_z,
                                                    d_m_joint_coordinate_x,
                                                    d_m_joint_coordinate_y,
                                                    d_m_joint_coordinate_z,

                                                    // STATISTICAL Parameters
                                                    d_m_is_deposited_particle,
                                                    d_m_is_escaped_particle,
                                                    d_m_is_error_particle,
                                                    d_m_is_stagnant_particle,
                                                    d_m_is_ghost_particle,
                                                    N_Cells,
                                                    N_DOF,
                                                    N_Particles_released,
                                                    time_step);

    cudaDeviceSynchronize();

    cudaEventCreate(&stop);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime, start, stop);
    cout << "[INFORMATION] Time taken for the kernel to execute : " << elapsedTime << " ms" << endl;

    // record the time taken for the transfer of data from device to host
    cudaEvent_t start1, stop1;
    cudaEventCreate(&start1);
    cudaEventRecord(start1, 0);


    // Lets transfer the position values back to the host
    checkCudaErrors(cudaMemcpy(position_X.data(), d_m_particle_position_x, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(position_Y.data(), d_m_particle_position_y, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(position_Z.data(), d_m_particle_position_z, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));

    // Lets transfer the previous position values back to the host
    checkCudaErrors(cudaMemcpy(position_X_old.data(), d_m_particle_previous_position_x, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(position_Y_old.data(), d_m_particle_previous_position_y, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(position_Z_old.data(), d_m_particle_previous_position_z, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));

    // Lets transfer the velocity values back to the host
    checkCudaErrors(cudaMemcpy(velocityX.data(), d_m_particle_velocity_x, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(velocityY.data(), d_m_particle_velocity_y, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(velocityZ.data(), d_m_particle_velocity_z, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));

    // transfer previous velocity values back to the host
    checkCudaErrors(cudaMemcpy(velocityX_old.data(), d_m_particle_previous_velocity_x, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(velocityY_old.data(), d_m_particle_previous_velocity_y, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(velocityZ_old.data(), d_m_particle_previous_velocity_z, N_Particles * sizeof(double), cudaMemcpyDeviceToHost));

    // Transfer the current cell and previous cell back to the host
    checkCudaErrors(cudaMemcpy(currentCell.data(), d_m_current_cell, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(previousCell.data(), d_m_previous_cell, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));

    // Transfer Statistical variables back to the host
    checkCudaErrors(cudaMemcpy(isParticleDeposited.data(), d_m_is_deposited_particle, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(isEscapedParticle.data(), d_m_is_escaped_particle, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(isStagnantParticle.data(), d_m_is_stagnant_particle, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(isErrorParticle.data(), d_m_is_error_particle, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(isGhostParticle.data(), d_m_is_ghost_particle, N_Particles * sizeof(int), cudaMemcpyDeviceToHost));

    // record the time taken for the transfer of data from device to host
    cudaEventCreate(&stop1);
    cudaEventRecord(stop1, 0);
    cudaEventSynchronize(stop1);

    float elapsedTime1;
    cudaEventElapsedTime(&elapsedTime1, start1, stop1);
    cout << "[INFORMATION] Time taken for the data transfer from device to host : " << elapsedTime1 << " ms" << endl;

}




