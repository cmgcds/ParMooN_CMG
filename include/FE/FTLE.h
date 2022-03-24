
// Idenitify the vertex neibhours
// 
#ifndef __FTLE__
#define __FTLE__

#include<vector>
#include<unordered_set>

class FTLE
{
    public:
        double* x_pos;  // Stores x position
        double* y_pos;  // Stores y position
        int* cellId;    // Stores the cell ID of all the cells;
        int N_Particles;
        int searchDepth;  // to how many levels the particles have to be tracked. 

        std::vector< std::unordered_set<int> > CellNeibhours;


        TFEVectFunct2D* VelocityVectFunction;

        // COnstructor
        FTLE(TFEVectFunct2D* FEVectFunction,int searchDepth);

        // Get the Velocity of the particle based on current position. 
        double* obtainVelocity();


};

#endif

// Need Each DOF to be mapped to every cell possible 
// DOF - 1 --> 
// Global DOF has all the cell mapping 


