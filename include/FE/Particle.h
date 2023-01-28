#include<vector>


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
        std::vector<int> cellsOnBoundary;
        std::vector<int> cellsOnBoundaryFaces;

        
        // Particle Co-ordinates 
        std::vector<double> position_X;
        std::vector<double> position_Y;
        std::vector<double> position_Z;

        // Particle Old
        std::vector<double> position_X_old;
        std::vector<double> position_Y_old;
        std::vector<double> position_Z_old;

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

        // Position 

        //-- constructor--//
        TParticles(int N_Particles, double circle_x, double circle_y, double radius, TFESpace3D *fespace);

        // -- Methods to Initialise Particles -- //
        void  Initialiseparticles(int N_Particles,double circle_x, double circle_y, double radius,TFESpace3D* fespace);

        void interpolateNewVelocity(double timeStep,TFEVectFunct3D* VelocityFEVectFunction);

        void OutputFile(const char* filename);
};

