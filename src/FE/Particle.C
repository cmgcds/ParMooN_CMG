#include <AllClasses.h>
#include <Domain.h>
#include <Database.h>
#include <DiscreteForm3D.h>
#include <FEDatabase3D.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <Assemble3D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <DeformMesh3D.h>
#include <LinAlg.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <FEVectFunct3D.h>
#include <Constants.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <NodalFunctional3D.h>
#include <BaseCell.h>
#include <Collection.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include <FTLE.h>
#include <stdlib.h>
#include <mkl.h>
#include <cmath>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <colors.h>

#include <omp.h>
#include <random>
#include<fstream>

#include <Particle.h>
#include<random>
#include <BoundFace.h>
#include <IsoJointEqN.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <FEFunction3D.h>
#include <InterfaceJoint3D.h>

/*
Jan-4 , Made changes to the particle constructor class to give which plane the inlet lies in
*/

struct VectorStruct
{
    double x;
    double y;
    double z;
};

/* Constructor for Particle Class */
TParticles::TParticles(int N_Particles_, double circle_x, double circle_y, double radius, TFESpace3D *fespace)
{
    // Initialise particle Stream
    this->N_Particles = N_Particles_;

    // Resize all the vectors
    DiameterParticles.resize(N_Particles, 10e-6);
    position_X.resize(N_Particles, 0.0);
    position_Y.resize(N_Particles, 0.0);
    position_Z.resize(N_Particles, 0.0);
    position_X_old.resize(N_Particles, 0.0);
    position_Y_old.resize(N_Particles, 0.0);
    position_Z_old.resize(N_Particles, 0.0);
    Density.resize(N_Particles, 1266);

    currentCell.resize(N_Particles, 0);
    previousCell.resize(N_Particles, 0);

    velocityX.resize(N_Particles, 0);
    velocityY.resize(N_Particles, 0);
    velocityZ.resize(N_Particles, 0);

    velocityX_old.resize(N_Particles, 0);
    velocityY_old.resize(N_Particles, 0);
    velocityZ_old.resize(N_Particles, 0);

    isParticleDeposited.resize(N_Particles,false);

    // Create N particles
    Initialiseparticles(N_Particles_, circle_x, circle_y, radius, fespace);
}

void TParticles::Initialiseparticles(int N_Particles, double circle_x, double circle_y, double radius, TFESpace3D *fespace)
{
    // Distribute particles within uniform Distribution

    double StartX = circle_x - radius;
    double EndX = circle_x + radius;
    double StartY = circle_y - radius;
    double EndY = circle_y + radius;

    int N_Cells = fespace->GetCollection()->GetN_Cells();

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(circle_x - radius,circle_y + radius);


    // For the siminhale, the inlet in in 

    for (int particleNo = 0; particleNo < N_Particles; particleNo++)
    {
        // Use rejection Sampling for points inside circle
        double y =  distribution(generator);
        double z =  distribution(generator);

        double distance = sqrt(y*y + z*z) ;

        int counter  = 0;
        //0.05 is added as saftey parameter to not to drop particles on surfaces. 
        while( distance > 0.95*radius  )
        {
            y =  distribution(generator);
            z =  distribution(generator);
            distance = sqrt(y*y + z*z) ;
            counter ++;

            if(counter % 10 == 0)
            {
                cout << "y : " << y << " z : " << z << " Dist : " << distance << " Counter : " << counter << " \n";
            }

        }
         
        // double random = radius * sqrt(rand());
        // double theta = rand() * 2 * Pi;

        // // Cartesian Co-ordinates
        // double z = circle_x + radius * cos(theta);
        // double y = circle_y + radius * sin(theta);

        position_X[particleNo] = y;
        position_Y[particleNo] = 0.001;
        position_Z[particleNo] = z;

        // Identify, which cell the particle Belongs to.
        int N_Cells = fespace->GetCollection()->GetN_Cells();
        for (int cellId = 0; cellId < N_Cells; cellId++)
        {
            TBaseCell *cell = fespace->GetCollection()->GetCell(cellId);
            bool insideDomain = false;
            if (cell->PointInCell(position_X[particleNo], position_Y[particleNo], position_Z[particleNo]))
            {
                insideDomain = true;
                currentCell[particleNo] = cellId;
                break;
            }
        }
        
    }

        
    // Find all the cells with atleast one boundary face on the Given Surface
    // Also tag the Joint ID's which belong to the Boundary edge 
    // Note : This approach has edge cases on the cells which might have two boundary Surfaces, 
    //        Which might result in picking either one of  the boundary surface on the cell.
    int MaxLen;
    int N_Joints;
    const int *TmpLen;
    const int *TmpFV;
    for (int cellNo = 0; cellNo < N_Cells; cellNo++)
    {
        // Get the cell
        TBaseCell *cell = fespace->GetCollection()->GetCell(cellNo);

        // NOw Obtain the faces of the cell.
        FE3D elementId = fespace->GetFE3D(cellNo, cell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        TBoundFace* Bdface;					// Pointer to Boundary Face in 3D Cell
        TBoundComp *BoundComp;
        N_Joints = cell->GetN_Joints();
        bool BoundaryJointFound = false;
        for (int jointId = 0; jointId < N_Joints; jointId++)
        {
            TJoint *Joint = cell->GetJoint(jointId);
            
            // Pick the cells with Boundary Surfaces and Ignore Inlet Boundary faces. 
            if (Joint->GetType() == BoundaryFace )
            {
                Bdface = (TBoundFace*)Joint;
				BoundComp = Bdface->GetBoundComp();
				int bdid = BoundComp->GetID();
                if(bdid != 0)
                {
                    cellsOnBoundary.push_back(cellNo);
                    cellsOnBoundaryFaces.push_back(jointId);
                    BoundaryJointFound = true;
                }
                
                
            }
            if(BoundaryJointFound) break;
        }
    }

    /// Print the Statistics of the Particle INitialisation 
    cout << " Number of particles Initialised : " << N_Particles <<endl;
    cout << " No ofBoundary cells Identified  : " << cellsOnBoundary.size()  <<endl;
    cout << " No ofBoundary faces Identified : " << cellsOnBoundaryFaces.size()  <<endl;
    
}

// Generates output file for Visualisation
void TParticles::OutputFile(const char* filename)
// Generates an CSV file of all co-ordinates for visualisation purposes. 
{
    // Obtain all the vector data and then write them as CSV files for visualisation
    std::ofstream file(filename);

    if(!file.is_open())
    {
        cout << " Error in Openning the file " << filename <<endl;
        cout << " File : Particle.C :: Function : OutputFile()" <<endl;
        cout << " exiting Process()" <<endl;
        exit(0);
    }



    file <<  "x"<<","<<"y" <<","<<"z"<<","<<"val"<<"\n";
    for ( int i = 0 ; i < position_X.size(); i++)
    {
        file  <<  position_X[i] <<"," << position_Y[i]<<","<< position_Z[i] <<","<<1<<"\n";
    }

    file.close();
}


// Helper Function

// Function to Interpolate the Velocity of the Particle in next time Step 
void TParticles::interpolateNewVelocity(double timeStep,TFEVectFunct3D* VelocityFEVectFunction)
{
    // Constants required for computation 
    double densityFluid = 1.1839; 
    double densityParticle= 1266;
    double g_x = 0;
    double g_y = 0;
    // double g_z = 1.0 / 51.6414;
    double g_z = 9.81;

    double dynamicViscosityFluid = 0.0000186;
    double lambda  = 0.00000007;

    // Particle Diameter = 8 Micrometers
    double particleDiameter = 2.5e-6;
    
    double mass_particle    = (Pi*pow(particleDiameter,3)*densityParticle)/6;
    double mass_div_dia = mass_particle/particleDiameter;
    // For the First Term
    double intertialConstant = (3./4.)*(densityFluid/densityParticle)*(1/particleDiameter);


    // For the second term
    double gForceConst_x = g_x;
    double gForceConst_y = g_y;
    double gForceConst_z = g_z;

    int MaxLen;
    int N_Joints;
    const int *TmpLen;
    const int *TmpFV;

    int ErrorParticles = 0;

    int FirstTime =1;

    //Lambda Function to compute CD/CC
    auto CD_CC = [&](double particleVel, double fluidVel)
    {
        double Re_Particle = densityFluid * particleDiameter* fabs(fluidVel - particleVel) / dynamicViscosityFluid;
        double CD = (24/Re_Particle) * ( 1 + 0.15*pow(Re_Particle,0.687));
        double CC =  1.0 + ((2*lambda)/particleDiameter) * ( 1.257 + 0.4*exp(-1.0 * ((1.1 * particleDiameter)/(2*lambda) ) ) );
        CC =1.0;
        // return CD/CC;
        return 1;

    };

    for ( int i = 0 ; i < N_Particles; i++)
    {   
        // cout << " =====================================================================================================================  " <<endl;
        if (isParticleDeposited[i] == true) continue;
        double values[4];
        int CellNo = currentCell[i];
        TFESpace3D* fespace = VelocityFEVectFunction->GetComponent(0)->GetFESpace3D();
        TBaseCell* cell = fespace->GetCollection()->GetCell(CellNo);
        VelocityFEVectFunction->GetComponent(0)->FindGradientLocal(cell,CellNo,position_X[i],position_Y[i],position_Z[i],values);
        double fluidVelocityX = values[0];
        VelocityFEVectFunction->GetComponent(1)->FindGradientLocal(cell,CellNo,position_X[i],position_Y[i],position_Z[i],values);
        double fluidVelocityY = values[0];
        VelocityFEVectFunction->GetComponent(2)->FindGradientLocal(cell,CellNo,position_X[i],position_Y[i],position_Z[i],values);
        double fluidVelocityZ = values[0];
        
        double cdcc_x = CD_CC(velocityX[i],fluidVelocityX);
        double cdcc_y = CD_CC(velocityY[i],fluidVelocityY);
        double cdcc_z = CD_CC(velocityZ[i],fluidVelocityZ);

        // The RHS will be 
        double rhs_x  = intertialConstant *cdcc_x* fabs(fluidVelocityX - velocityX[i])*(fluidVelocityX - velocityX[i]) + gForceConst_x*(densityFluid - densityParticle)/densityParticle;
        double rhs_y  = intertialConstant *cdcc_y* fabs(fluidVelocityY - velocityY[i])*(fluidVelocityY - velocityY[i]) + gForceConst_y*(densityFluid - densityParticle)/densityParticle;
        double rhs_z  = intertialConstant *cdcc_z* fabs(fluidVelocityZ - velocityZ[i])*(fluidVelocityZ - velocityZ[i]) + gForceConst_z*(densityFluid - densityParticle)/densityParticle;

        // cout << "-- intertialConstant : " << intertialConstant   << "  CDCC : " << (CD_CC(velocityX[i],fluidVelocityX)) << " fluidVelocityX: " << fluidVelocityX << "    velocityX : " <<velocityX[i] << "\n";
        // cout << "-- intertialConstant     : " << intertialConstant <<"  CDCC : " << (CD_CC(velocityY[i],fluidVelocityY)) <<"     fluidVelocityY: " << fluidVelocityY << "    velocityY : " <<velocityY[i] << "\n";
        // cout << "-- intertialConstant : " << intertialConstant << "  CDCC : " << (CD_CC(velocityZ[i],fluidVelocityZ)) <<" fluidVelocityZ: " << fluidVelocityZ << "    velocityZ : " <<velocityZ[i] << "\n";
        
        // Now, Compute the Updated Velocity of particle Using Forward Euler

        double velocityParticle_X_New = rhs_x * (timeStep) + velocityX[i];
        double velocityParticle_Y_New = rhs_y * (timeStep) + velocityY[i];
        double velocityParticle_Z_New = rhs_z * (timeStep) + velocityZ[i];

        // cout << "-- vel new x: " << velocityParticle_X_New << " rhs: " << rhs_x << " t: " << (timeStep) << "\n";
        // cout << "-- vel new y: " << velocityParticle_Y_New << " rhs: " << rhs_y << " t : " << (timeStep) << "\n";
        // cout << "-- vel new z: " << velocityParticle_Z_New << " rhs: " << rhs_z << " t : " << (timeStep) << "\n";
        
        // Transfer current velocity as velocity old
        velocityX_old[i] = velocityX[i];
        velocityY_old[i] = velocityY[i];
        velocityZ_old[i] = velocityZ[i];
        
        // Udpate new velocity as current velocity
        velocityX[i] = velocityParticle_X_New;
        velocityY[i] = velocityParticle_Y_New;
        velocityZ[i] = velocityParticle_Z_New;

        // cout << "-- " << "Vel x : " << velocityX[i] <<  " old vel x :" << velocityX_old[i] << " \n";
        // cout << "-- " << "Vel y : " << velocityY[i] <<  " old vel y :" << velocityY_old[i] << " \n";
        // cout << "-- " << "Vel z : " << velocityZ[i] <<  " old vel z :" << velocityZ_old[i] << " \n";

        //  Now Compute the new position of the particle using the updated velocity
        position_X_old[i] = position_X[i];
        position_Y_old[i] = position_Y[i];
        position_Z_old[i] = position_Z[i];

       

        position_X[i] += timeStep * 0.5 * (velocityX_old[i] + velocityX[i]);
        position_Y[i] += timeStep * 0.5 * (velocityY_old[i] + velocityY[i]);
        position_Z[i] += timeStep * 0.5 * (velocityZ_old[i] + velocityZ[i]);

        // cout <<" Old Position x : " << position_X_old[i] << " New Position x: " << position_X[i] <<endl;
        // cout <<" Old Position y: " << position_Y_old[i] << " New Position y: " << position_Y[i] <<endl;
        // cout <<" Old Position z: " << position_Z_old[i] << " New Position z: " << position_Z[i] <<endl;

        // Update the current position of cell. 
        // Check if the particle Exists in the current Domain
        int N_Cells = fespace->GetN_Cells();
        bool insideDomain = false;
        for (int cellId = 0; cellId <  N_Cells; cellId++)
        {
            cell = fespace->GetCollection()->GetCell(cellId);

            if (cell->PointInCell(position_X[i], position_Y[i],position_Z[i]))
            {
                insideDomain = true;
                previousCell[i] = currentCell[i];
                currentCell[i] = cellId;
                break;
            }
        }

        // If Not, particle is deposited. 
        if(!insideDomain)
        {
            // Find the last Cell of the particle 
            int cellNo = currentCell[i];
            // int  jointID= cellsOnBoundaryFaces[i];
            
            TBaseCell* cell = fespace->GetCollection()->GetCell(cellNo);

            FE3D elementId = fespace->GetFE3D(cellNo, cell);
            TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
            TFEDesc3D *fedesc = element->GetFEDesc3D();
            cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
            int** JointDOFs = fedesc->GetJointDOF();
            

            TBoundFace* Bdface;					// Pointer to Boundary Face in 3D Cell
            TBoundComp *BoundComp;
            int N_Joints = cell->GetN_Joints();
            
            // Identify the Boundary Face
            int jointID = -999999999;
            for (int jj = 0; jj < N_Joints; jj++)
            {
                TJoint *Joint = cell->GetJoint(jj);
                
                // Pick the cells with Boundary Surfaces and Ignore Inlet Boundary faces. 
                if (Joint->GetType() == BoundaryFace )
                {
                    Bdface = (TBoundFace*)Joint;
                    BoundComp = Bdface->GetBoundComp();
                    int bdid = BoundComp->GetID();
                    if(bdid != 0)
                    {
                       jointID = jj;
                       break;
                    }
                    
                    
                }
            }
 

            if(jointID == -999999999) // Last cell was not a Boundary cell
            {
                // Check for Boundary of the neighbours for particles. 
                cout << "Error" <<endl;
                isParticleDeposited[i] = true;
                m_ErrorParticlesCount++;
                break;

            }


            TJoint *Joint = cell->GetJoint(jointID); 
            double x1, x2, x3,y1,y2,y3,z1,z2,z3;
            
            cell->GetVertex(TmpFV[jointID*MaxLen+0])->GetCoords(x1, y1, z1);
            cell->GetVertex(TmpFV[jointID*MaxLen+1])->GetCoords(x2, y2, z2);
            double t11 = x2-x1; 
            double t12 = y2-y1; 
            double t13 = z2-z1;
            double len = sqrt(t11*t11 + t12*t12 + t13*t13);
            t11 /= len; t12 /= len; t13 /= len;

            cell->GetVertex(TmpFV[jointID*MaxLen+(TmpLen[jointID]-1)])->GetCoords(x2, y2, z2);
            double t21 = x2-x1; double t22 = y2-y1;double t23 = z2-z1;
            len = sqrt(t21*t21 + t22*t22 + t23*t23);
            t21 /= len; t22 /= len; t23 /= len;

            double N1 = t12*t23 - t13*t22;
            double N2 = t13*t21 - t11*t23;
            double N3 = t11*t22 - t12*t21;
            len = sqrt(N1*N1 + N2*N2 + N3*N3);
            N1 /= len; N2 /= len; N3 /= len;
            
            
            VectorStruct firstPoint;
            VectorStruct secondPoint;
            VectorStruct LineVector;
            VectorStruct pointOnSurface;
            VectorStruct normalSurface;
            VectorStruct temp1;
            VectorStruct temp2;


            firstPoint.x = position_X_old[i];firstPoint.y = position_Y_old[i];firstPoint.z = position_Z_old[i];
            secondPoint.x = position_X[i];secondPoint.y = position_Y[i];secondPoint.z = position_Z[i];

            normalSurface.x = N1;normalSurface.y = N2;normalSurface.z = N3;
            pointOnSurface.x = x1;pointOnSurface.y = y1;pointOnSurface.z = z1; 

            //  u = p1 - p0
            LineVector.x = firstPoint.x - secondPoint.x;
            LineVector.y = firstPoint.y - secondPoint.y;
            LineVector.z = firstPoint.z - secondPoint.z;
            

            // Dot 
            double dot =  normalSurface.x*firstPoint.x + normalSurface.y*firstPoint.y + normalSurface.z*firstPoint.z;

            if(fabs(dot - 0.0) > 1e-3)
            {
                // w = p0 - pC0
                temp1.x = firstPoint.x - pointOnSurface.x;
                temp1.y = firstPoint.y - pointOnSurface.y;
                temp1.z = firstPoint.z - pointOnSurface.z;

                double fac = -1.0* (normalSurface.x*temp1.x + normalSurface.y*temp1.y + normalSurface.z*temp1.z);
                fac/= fac;

                // u = u*fac
                LineVector.x = LineVector.x*fac;
                LineVector.y = LineVector.y*fac;
                LineVector.z = LineVector.z*fac;

                temp2.x = LineVector.x + firstPoint.x;
                temp2.y = LineVector.y + firstPoint.y;
                temp2.z = LineVector.z + firstPoint.z;

                // Mark the Particle as deposited 
                isParticleDeposited[i] = true;

                position_X[i] = temp2.x;
                position_Y[i] = temp2.y;
                position_Z[i] = temp2.z;
            }
            else
            {
                // Mark the Particle as deposited 
                isParticleDeposited[i] = true;
                // Ghost particle, Make the vertex as deposition 
                position_X[i] = x1;
                position_Y[i] = y1;
                position_Z[i] = z1;
            }
        }


        // check if the particle is in any border cell


        
        
    }   
    
    cout << " Difference in position X : " << Ddot(N_Particles,position_X.data(),position_X_old.data()) << "\n";
    cout << " Difference in position Y : " << Ddot(N_Particles,position_Y.data(),position_Y_old.data()) << "\n";
    cout << " Difference in position Z : " << Ddot(N_Particles,position_Z.data(),position_Z_old.data()) << "\n";

    

    int depositedCount = 0;
    for (int l = 0 ; l < isParticleDeposited.size(); l++)
        if(isParticleDeposited[l] == true)
            depositedCount++;
    int NotDeposited = N_Particles - depositedCount;
    cout << "No of particles Deposited : " << depositedCount <<endl;
    cout << "percentage of particles Deposited : " << ( double(N_Particles -  depositedCount)/(double)N_Particles ) * 100 << " % " <<endl;
    cout << "Error particles Accumulaated : " << m_ErrorParticlesCount <<endl;


}





