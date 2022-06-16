

// #include<vector>
#include <AllClasses.h>
#include<FEVectFunct2D.h>
#include<FESpace2D.h>
#include <vector>
#include<set>
#include <Output2D.h>

#ifndef __FTLE__
#define __FTLE__

class FTLE
{
    public:
        

        int N_Particles;
        int N_cells;
        int searchDepth;  // to how many levels the particles have to be tracked. 


        // Array Structures 
        double* x_pos;  // Stores x position
        double* y_pos;  // Stores y position

        double* x_pos_initial;          // Stores x position
        double* y_pos_initial;          // Stores y position
        int* cellIds;                   // Stores the cell ID of all the cells;
        int* isParticleOnBoundary;     // Stores if the particle is present on the boundary or not. 

        
        std::vector< std::vector<int> > Cell2DOFMapping;
        std::vector< std::vector<int> > DOF2CellMapping;

        std::vector< std::vector<int> > DOF_Neibhours;
        std::vector< std::vector<int> > DOF_Neibhours_depth;
        std::vector<bool> particleInsideDomain;

        std::vector<int> neibhourLeft;
        std::vector<int> neibhourRight;
        std::vector<int> neibhourTop;
        std::vector<int> neibhourBottom;


        TFEVectFunct2D* VelocityVectFunction;
        TFESpace2D* fespace;
        TFEVectFunct2D* u_Vectfunc_Combined;
        TFEVectFunct2D* v_Vectfunc_Combined;


        std::vector<int> onBoundary;

        std::vector<double> FTLEValues;

        
        

        // COnstructor
        FTLE(TFESpace2D* ftleFespace,TFEVectFunct2D* FEVectFunction,int searchDepth,TFEVectFunct2D* u_Vectfunc_Combined,TFEVectFunct2D* v_Vectfunc_Combined );

        // Get the Velocity of the particle based on current position. 
        double* computeFTLE(double timeStep, int T, int startT);

        // Write the FTLE valyes to a VTK, So that they can be visualised using paraview
        void write_vtk(std::string fileName);


};

#endif

// Need Each DOF to be mapped to every cell possible 
// DOF - 1 --> 
// Global DOF has all the cell mapping 


