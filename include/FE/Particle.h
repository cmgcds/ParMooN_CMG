#include<vector>
#include<string>
#include <map>
#include <vector>
#include <AllClasses.h>
#include <FEFunction3D.h>

#ifndef __Particle__
#define __Particle__

// Class of Mono Disperse Particles
class TParticles
{
    public:
        // int Num of Particles 
        int N_Particles;

        // Size(diameter) of particles
        std::vector<double> m_particle_diameter;

        // Density of the Moecules 
        std::vector<double> m_particle_density;

        // Gravity
        double m_gravity_x;
        double m_gravity_y;
        double m_gravity_z;

        // mean free path 
        double m_lambda;

        // fluid density
        double m_fluid_density;

        // Dynamic Viscosity
        double m_fluid_dynamic_viscosity;


        // mass of the particles 
        std::vector<double> MassParticles;

        // Bool for deposition
        std::vector<int> isParticleDeposited;

        //  Track cells with Boundary Faces. 
        // Generate a map of cells with boundary faces using map 
        std::map<int,std::vector<int> > m_mapBoundaryFaceIds;
        std::map<int,int > m_cornerTypeOfBoundCells;
        std::map<int,int > m_jointidOfBoundCells;

        
        std::map<int,std::vector<double> > m_BoundaryDOFsOnCell; 
        // Track the Neibhour DOF's with a vector
        


        // Particle Co-ordinates 
        std::vector<double> position_X;
        std::vector<double> position_Y;
        std::vector<double> position_Z;

        // Particle Old
        std::vector<double> position_X_old;
        std::vector<double> position_Y_old;
        std::vector<double> position_Z_old;

        // Previous Particle Co-ordinates for stagnant check
        std::vector<double> previousPosition_X;
        std::vector<double> previousPosition_Y;
        std::vector<double> previousPosition_Z;

        //Particle Cell information
        std::vector<int> currentCell;
        std::vector<int> previousCell;

        // Particle velocity Information 
        std::vector<double> velocityX_old;
        std::vector<double> velocityY_old;
        std::vector<double> velocityZ_old;

        std::vector<double> velocityX;
        std::vector<double> velocityY;
        std::vector<double> velocityZ;


        //Error Particles Count
        int m_ErrorParticlesCount = 0;
        std::vector<int> isErrorParticle;

        // Escaped Particles Count
        int m_EscapedParticlesCount = 0;
        std::vector<int> isEscapedParticle;

				// Stagnant Particles Count
				int m_StagnantParticlesCount = 0;
				std::vector<int> isStagnantParticle;


        // Number of particles released
        int m_ParticlesReleased = 0;

        // Number of Ghose Particles 
        int m_ghostParticlesCount = 0;
        std::vector<int>  isGhostParticle;
        // Position 

        //-- constructor--//
	TParticles(int N_Particles, double circle_x, double circle_y, double radius, TFESpace3D *fespace, std::string resumeFile, int *StartNo);

        // -- Methods to Initialise Particles -- //
        void  Initialiseparticles(int N_Particles,double circle_x, double circle_y, double radius,TFESpace3D* fespace);

        // -- Methods to Update Simulation Parmameters -- //
        void InitialiseParticleParameters(int N_Particles_);

        void interpolateNewVelocity(double timeStep,TFEVectFunct3D* VelocityFEVectFunction, TFESpace3D* fespace);

        void interpolateNewVelocity_Parallel(double timeStep,TFEVectFunct3D* VelocityFEVectFunction, TFESpace3D* fespace);

        void OutputFile(const char* filename);

				int UpdateParticleDetailsFromFile(std::string filename);

				void printUpdatedParticleDetailStats();

				void detectStagnantParticles();

        #ifdef _CUDA

          // FEM Data Structures
          // Cell Vertex  Co-ordinates
          double* h_m_cell_vertices_x;
          double* h_m_cell_vertices_y;
          double* h_m_cell_vertices_z;
          // DOF Indices
          int* h_m_global_dof_indices;
          int* h_m_begin_indices;


          // velocity Arrays
          double* h_m_velocityX;
          double* h_m_velocityY;
          double* h_m_velocityZ;

          // n_basis_functions for each cell
          int h_m_n_basis_functions;

          // basis functions for each cell
          double* h_m_basis_functions_values;

          // To store the timestep 
          // NOte : This value is not used in actual host rotuine, rather these value are 
          // directly passed as fnction parameter
          double h_m_time_step;


          // -- Depositon Related Datastructures -- //

          // To save if a cell is boundary cell or not
          int* h_m_is_boundary_cell;

          // To save the corner id of the given cell
          int* h_m_corner_id;

          // to check if the cell has boundary DOF or not
          int* h_m_is_boundary_dof_present;

          // to Store the x-coordinates of the boundary DOF's
          double* h_m_boundary_dof_x;

          // to Store the y-coordinates of the boundary DOF's
          double* h_m_boundary_dof_y;

          // to Store the z-coordinates of the boundary DOF's
          double* h_m_boundary_dof_z;

          // Store the Joint id of the boundary cell
          int* h_m_joint_id;

          // Store the normal in x direction of the boundary cell
          double* h_m_joint_normal_x;

          // Store the normal in y direction of the boundary cell
          double* h_m_joint_normal_y;

          // Store the normal in z direction of the boundary cell
          double* h_m_joint_normal_z;

          // Store the Co-ordinate of the boundary joint - x direction
          double* h_m_joint_coordinate_x;

          // Store the Co-ordinate of the boundary joint - y direction
          double* h_m_joint_coordinate_y;

          // Store the Co-ordinate of the boundary joint - z direction
          double* h_m_joint_coordinate_z;


          // ------------ DEVICE VARIABLES DECLARATION -------------- //
          // -- FE Related Datastructures -- //
          // Lets delcare all the cuda variables 
          // Cell Co-ordinates
          double* d_m_cell_vertices_x;
          double* d_m_cell_vertices_y;
          double* d_m_cell_vertices_z;

          // DOF Indices
          int* d_m_global_dof_indices;
          int* d_m_begin_indices;


          // velocity Arrays
          double* d_m_velocity_nodal_values_x;
          double* d_m_velocity_nodal_values_y;
          double* d_m_velocity_nodal_values_z;

          // copy the current cell location to GPU
          int* d_m_current_cell;

          // copy the previous cell location to GPU
          int* d_m_previous_cell;

          // Allocate Arrays for the particle position
          double* d_m_particle_position_x;
          double* d_m_particle_position_y;
          double* d_m_particle_position_z;

          // Allocate Arrays for the particle previous position
          double* d_m_particle_previous_position_x;
          double* d_m_particle_previous_position_y;
          double* d_m_particle_previous_position_z;

          // Allocate Arrays for the particle previous position for stagnant particles
          double* d_m_particle_stagnant_position_x;
          double* d_m_particle_stagnant_position_y;
          double* d_m_particle_stagnant_position_z;

          // Allocate memory for the particle velocity
          double* d_m_particle_velocity_x;
          double* d_m_particle_velocity_y;
          double* d_m_particle_velocity_z;

          // Allocate memory for the particle previous velocity
          double* d_m_particle_previous_velocity_x;
          double* d_m_particle_previous_velocity_y;
          double* d_m_particle_previous_velocity_z;

          // Allocate memory to store the interpolated velocity values at the current particle position
          double* d_m_interpolated_velocity_x;
          double* d_m_interpolated_velocity_y;
          double* d_m_interpolated_velocity_z;
          
          // Allocate memory for the density of particles
          double* d_m_particle_density;

          // Allocate memory for the diameter of particles
          double* d_m_particle_diameter;

          // Fluid density 
          double* d_m_fluid_density;

          // Fluid viscosity
          double* d_m_dynamic_viscosity_fluid;

          // lambda
          double* d_m_lambda;

          // Gravity parameters
          double* d_m_gravity_x;
          double* d_m_gravity_y;
          double* d_m_gravity_z;
          
          // Store the time step
          double* d_m_time_step;


          // -- Deposition Related Datastructures -- //
          // To save if a cell is boundary cell or not
          int* d_m_is_boundary_cell;

          // To save the corner id of the given cell
          int* d_m_corner_id;

          // to check if the cell has boundary DOF or not
          int* d_m_is_boundary_dof_present;

          // to Store the x-coordinates of the boundary DOF's
          double* d_m_boundary_dof_x;

          // to Store the y-coordinates of the boundary DOF's
          double* d_m_boundary_dof_y;

          // to Store the z-coordinates of the boundary DOF's
          double* d_m_boundary_dof_z;

          // Store the Joint id of the boundary cell
          int* d_m_joint_id;

          // Store the normal in x direction of the boundary cell
          double* d_m_joint_normal_x;

          // Store the normal in y direction of the boundary cell
          double* d_m_joint_normal_y;

          // Store the normal in z direction of the boundary cell
          double* d_m_joint_normal_z;

          // Store the Co-ordinate of the boundary joint - x direction
          double* d_m_joint_coordinate_x;

          // Store the Co-ordinate of the boundary joint - y direction
          double* d_m_joint_coordinate_y;

          // Store the Co-ordinate of the boundary joint - z direction
          double* d_m_joint_coordinate_z;


          // To save the Summary Statistics of a given particle. 
          
          // escaped particles status
          int* d_m_is_escaped_particle;

          // error particles status
          int* d_m_is_error_particle;

          // stagnant particles status
          int* d_m_is_stagnant_particle;

          // ghost particles status
          int* d_m_is_ghost_particle;

          // deposited particles status
          int* d_m_is_deposited_particle;


          // --- CUDA RELATED FUNCTIONS -- //
          void  CD_CC_Cuda();
          void  FindValueLocal_Parallel();

          // Setup the essential data structures for CUDA
          void SetupCudaDataStructures(TFESpace3D* fespace);

          // Wrapper function for transfering velocity and N_particles_released to GPU
          void SetupVelocityValues(double* velocity_x_data, double* velocity_y_data, double* velocity_z_data, int N_Particles_released, int n_dof);

          // Host wrapper for calling interpolate function
          void InterpolateVelocityHostWrapper(double timeStep,int N_Particles_released, int N_DOF, int N_cells);

          // Host wrapper for detecting stagnant particles
          void DetectStagnantParticlesHostWrapper(int N_Particles_released);

        #endif


		protected:
				bool isStagnant(int i);
};

#endif
