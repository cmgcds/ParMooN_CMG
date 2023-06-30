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
        std::vector<double> DiameterParticles;

        // Density of the Moecules 
        std::vector<double> Density;

        // mass of the particles 
        std::vector<double> MassParticles;

        // Bool for deposition
        std::vector<bool> isParticleDeposited;

        //  Track cells with Boundary Faces. 
        // Generate a map of cells with boundary faces using map 
        std::map<int,std::vector<int> > m_mapBoundaryFaceIds;
        std::map<int,int > m_cornerTypeOfBoundCells;
        
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
        // Position 

        //-- constructor--//
        TParticles(int N_Particles, double circle_x, double circle_y, double radius, TFESpace3D *fespace);

        // -- Methods to Initialise Particles -- //
        void  Initialiseparticles(int N_Particles,double circle_x, double circle_y, double radius,TFESpace3D* fespace);

        void interpolateNewVelocity(double timeStep,TFEVectFunct3D* VelocityFEVectFunction, TFESpace3D* fespace);

        void interpolateNewVelocity_Parallel(double timeStep,TFEVectFunct3D* VelocityFEVectFunction, TFESpace3D* fespace);

        void OutputFile(const char* filename);

				int UpdateParticleDetailsFromFile(std::string filename);

				void printUpdatedParticleDetailStats();

				void detectStagnantParticles();

        #ifdef _CUDA
        
          void  CD_CC_Cuda();
          void  FindValueLocal_Parallel();

          // Setup the essential data structures for CUDA
          void SetupCudaDataStructures(TFESpace3D* fespace);
        
          // Particle Co-ordinates
          double* h_m_cell_vertices;
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

          // Lets delcare all the cuda variables 
          // Particle Co-ordinates
          double* c_m_cell_vertices;
          // DOF Indices
          int* c_m_global_dof_indices;
          int* c_m_begin_indices;

          // velocity Arrays
          double* c_m_velocityX;
          double* c_m_velocityY;
          double* c_m_velocityZ;

          // n_basis_functions for each cell
          int* c_m_n_basis_functions;

          // basis functions for each cell
          double* c_m_basis_functions_values;

        // 

        #endif


		protected:
				bool isStagnant(int i);
};

#endif