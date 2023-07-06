#include<vector>
#include<string>

// Class of Mono Disperse Particles
class TParticles
{
    public:
        // int Num of Particles 
        int N_Particles;

        // Size(diameter) of particles
        std::vector<double> DiameterParticles;

        // Density of the Molecules 
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

		protected:
				bool isStagnant(int i);
};
