

// #include<vector>
#include <AllClasses.h>
#include<FEVectFunct2D.h>
#include<FESpace2D.h>
#include <vector>
#include<set>

#ifndef __FTLE__
#define __FTLE__

class FTLE
{
    public:
        double* x_pos;  // Stores x position
        double* y_pos;  // Stores y position

        double* x_pos_initial;  // Stores x position
        double* y_pos_initial;  // Stores y position

        int* cellIds;    // Stores the cell ID of all the cells;
        int N_Particles;
        int searchDepth;  // to how many levels the particles have to be tracked. 
        
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

        

        // COnstructor
        FTLE(TFEVectFunct2D* FEVectFunction,int searchDepth,TFEVectFunct2D* uVect_combined,TFEVectFunct2D* vVect_combined );

        // Get the Velocity of the particle based on current position. 
        double* computeFTLE(double* uVelocityatPoint, double* vVelocityatPoint,double timeStep,double T,double startTime);




};

#endif

// Need Each DOF to be mapped to every cell possible 
// DOF - 1 --> 
// Global DOF has all the cell mapping 


