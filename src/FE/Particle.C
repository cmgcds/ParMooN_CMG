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
#include <string.h>
#include <sstream>


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

// #include <FTLE.h>
#include <stdlib.h>
#include <mkl.h>
#include <cmath>
#include <set>
#include <unordered_set>
#include <algorithm>
// #include <colors.h>

#include <omp.h>
#include <random>
#include <fstream>

#include <Particle.h>
#include <random>
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
    DiameterParticles.resize(N_Particles, 4.3e-6); // earlier 10e-6
    position_X.resize(N_Particles, 0.0);
    position_Y.resize(N_Particles, 0.0);
    position_Z.resize(N_Particles, 0.0);
    position_X_old.resize(N_Particles, 0.0);
    position_Y_old.resize(N_Particles, 0.0);
    position_Z_old.resize(N_Particles, 0.0);
    previousPosition_X.resize(N_Particles, 0.0);
    previousPosition_Y.resize(N_Particles, 0.0);
    previousPosition_Z.resize(N_Particles, 0.0);
    isErrorParticle.resize(N_Particles, 0);
    isEscapedParticle.resize(N_Particles, 0);
    isStagnantParticle.resize(N_Particles, 0);
    Density.resize(N_Particles, 914); // earlier 1266

		// initialise previousPositions to zero
		for (int i = 0; i < N_Particles; i++) {
			previousPosition_X[i] = 0.0;
			previousPosition_Y[i] = 0.0;
			previousPosition_Z[i] = 0.0;
		}

    currentCell.resize(N_Particles, 0);
    previousCell.resize(N_Particles, 0);

    velocityX.resize(N_Particles, 0);
    velocityY.resize(N_Particles, 0);
    velocityZ.resize(N_Particles, 0);

    velocityX_old.resize(N_Particles, 0);
    velocityY_old.resize(N_Particles, 0);
    velocityZ_old.resize(N_Particles, 0);

    isParticleDeposited.resize(N_Particles, false);

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
    std::uniform_real_distribution<double> distribution(circle_x - radius, circle_y + radius);

    // For the siminhale, the inlet in in k0

    for (int particleNo = 0; particleNo < N_Particles; particleNo++)
    {
        if(particleNo % 1000 == 0)
            cout << "Particle : " << particleNo <<endl;
        // Use rejection Sampling for points inside circle
        double y = distribution(generator);
        double z = distribution(generator);

        double distance = sqrt(y * y + z * z);

        int counter = 0;
        // 0.05 is added as saftey parameter to not to drop particles on surfaces.
        while (distance > 0.95 * radius)
        {
            y = distribution(generator);
            z = distribution(generator);
            distance = sqrt(y * y + z * z);
            counter++;

            if (counter % 10 == 0)
            {
                cout << "y : " << y << " z : " << z << " Dist : " << distance << " Counter : " << counter << " \n";
            }
        }



        position_X[particleNo] = y;
        position_Y[particleNo] = 0.001;
        position_Z[particleNo] = z;

        
       

        // Identify, which cell the particle Belongs to.
        int N_Cells = fespace->GetCollection()->GetN_Cells();
				int num_threads = (int) ceil(0.9 * omp_get_max_threads());
        #pragma omp parallel for num_threads(num_threads) schedule(static,1)  shared(N_Cells, fespace, particleNo,currentCell) 
        for (int cellId = 0; cellId < N_Cells; cellId++)
        {
            TBaseCell *cell = fespace->GetCollection()->GetCell(cellId);
            bool insideDomain = false;
            if (cell->PointInCell_Parallel(position_X[particleNo], position_Y[particleNo], position_Z[particleNo]))
            {
                insideDomain = true;
                currentCell[particleNo] = cellId;
                // cout <<" Particle : " << particleNo << " Inside Cell : " << cellId <<endl;
                #pragma omp cancel for

            }
        }
    }
     

    cout << " All Particles Initialised " << endl;


    // print the particle positions and their cell Ids using setwidth

    // for(int i = 0; i < N_Particles; i++)
    // {
    //     cout << "Particle : " << setw(3) << i << " x : " << setw(12) << position_X[i] << " y : " << setw(12) << position_Y[i] << " z : " << setw(12) << position_Z[i] << " Cell : " << setw(3) << currentCell[i] << endl;
    // }


    // ----- READ  PARTICLE ----------------- //

    // Find all the cells with atleast one boundary face on the Given Surface
    // Also tag the Joint ID's which belong to the Boundary edge
    // Note : This approach has edge cases on the cells which might have two boundary Surfaces,
    //        Which might result in picking either one of  the boundary surface on the cell.
    int MaxLen;
    int N_Joints;
    const int *TmpLen;
    const int *TmpFV;
    

    // Variables For Storing the Corner Parameters
    int N_corner_20 = 0;
    int N_corner_21 = 0;
    int N_corner_2 = 0;
    int N_corner_1 = 0;
    int N_corner_0 = 0;
    int N_corner_m1 = 0;
    // for (int cellNo = 0; cellNo < N_Cells; cellNo++)
		int num_threads = (int) ceil(0.9 * omp_get_max_threads());
    #pragma omp parallel for num_threads(num_threads) schedule(static,1)  shared(N_Cells, fespace,m_cornerTypeOfBoundCells, m_mapBoundaryFaceIds,m_BoundaryDOFsOnCell)
    for (int cellNo = 0; cellNo < N_Cells; cellNo++)
    {

        // Get the cell
        TBaseCell *cell = fespace->GetCollection()->GetCell(cellNo);

        // NOw Obtain the faces of the cell.
        FE3D elementId = fespace->GetFE3D(cellNo, cell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        // cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        TBoundFace *Bdface; // Pointer to Boundary Face in 3D Cell
        TBoundComp *BoundComp;
        N_Joints = cell->GetN_Joints();
        bool BoundaryJointFound = false;
        std::vector<int> tmp_face_ids;
        for (int jointId = 0; jointId < N_Joints; jointId++)
        {
            TJoint *Joint = cell->GetJoint(jointId);

            // Pick the cells with Boundary Surfaces and Ignore Inlet Boundary faces.
            if (Joint->GetType() == BoundaryFace)
            {
                Bdface = (TBoundFace *)Joint;
                BoundComp = Bdface->GetBoundComp();
                int bdid = BoundComp->GetID();
                // cellsOnBoundary.push_back(cellNo);
                tmp_face_ids.push_back(bdid);
                BoundaryJointFound = true;
            }
        }
        
        int cornerType = 0;

        
        // Assign the Face ids to the the map with cellNo as key

        // if tmp_face_ids has 0 and 2, its inlet corner and assign value of cornerType as 20
        // if tmp_face_ids has 1 and 2, its inlet corner and assign value of cornerType as 21
        // if tmp_face_ids have only 1, its inlet and assign value of cornerType as 1
        // if tmp_face_ids have only 0, its inlet and assign value of cornerType as 0
        // else its not inlet corner and assign value of cornerType as -1
        // Write the code using std::find
        if(BoundaryJointFound)
        {
            #pragma omp critical
            m_mapBoundaryFaceIds[cellNo] = tmp_face_ids;

            if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 0) != tmp_face_ids.end() && std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2) != tmp_face_ids.end())
            {
                cornerType = 20;
                N_corner_20++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 1) != tmp_face_ids.end() && std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2) != tmp_face_ids.end())
            {
                cornerType = 21;
                N_corner_21++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2) != tmp_face_ids.end())
            {
                cornerType = 2;
                N_corner_2++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 0) != tmp_face_ids.end())
            {
                cornerType = 0;
                N_corner_0++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 1) != tmp_face_ids.end())
            {
                cornerType = 1;
                N_corner_1++;
            }
            else
            {
                cornerType = -1;
                N_corner_m1++;
            }

            
            #pragma omp critical
            m_cornerTypeOfBoundCells[cellNo] = cornerType;
            

        }

        

        // In A given cell, Identify all the Boundary DOF Co-ordinates that exists within the cell. 
        // This is required to identify the deposition of the particle.
        std::vector<double> co_ordinates;
        int *globalDOFindex = fespace->GetGlobalNumbers();
        int *BeginIndex = fespace->GetBeginIndex();
        int start = BeginIndex[cellNo];
        int end = BeginIndex[cellNo + 1];
        int N_Active = fespace->GetActiveBound();
        
        for (int index = start; index < end; index++)
        {
            if (globalDOFindex[index] >= N_Active)
            {
                double xx = 0.;
                double yy = 0.;
                double zz = 0.;
                
                // mark the deposition as the vertex point.
                fespace->GetDOFPosition_Parallel(globalDOFindex[index], xx, yy, zz);
                co_ordinates.push_back(xx);
                co_ordinates.push_back(yy);
                co_ordinates.push_back(zz);

                // Assign the DOF Co-ordinates to the map with cellNo as key
                #pragma omp critical
                m_BoundaryDOFsOnCell[cellNo] = co_ordinates;

                 break;
            }
        }
    }

    cout <<" Reached here 2" <<endl;



    /// Print the Statistics of the Particle INitialisation
    cout << " Number of particles Initialised : " << N_Particles << endl;
    cout << " No ofBoundary cells Identified  : " << m_mapBoundaryFaceIds.size() << endl;

    // Print the N_count variables in a tablular form with fixed width with each parameter in a single line with all numbers right aligned
    std::cout << std::setw(10) << std::left << "N_corner_20" << " : " << std::setw(10) << std::right << N_corner_20 << std::endl;
    std::cout << std::setw(10) << std::left << "N_corner_21" << " : " << std::setw(10) << std::right << N_corner_21 << std::endl;
    std::cout << std::setw(10) << std::left << "N_corner_2" << " : " << std::setw(10) << std::right << N_corner_2 << std::endl;
    std::cout << std::setw(10) << std::left << "N_corner_1" << " : " << std::setw(10) << std::right << N_corner_1 << std::endl;
    std::cout << std::setw(10) << std::left << "N_corner_0" << " : " << std::setw(10) << std::right << N_corner_0 << std::endl;
    std::cout << std::setw(10) << std::left << "N_corner_m1" << " : " << std::setw(10) << std::right << N_corner_m1 << std::endl;
 

    // cout << " No ofBoundary faces Identified : " << Face_id_cellsOnBoundary.size() << endl;


}

// Generates output file for Visualisation
void TParticles::OutputFile(const char *filename)
// Generates an CSV file of all co-ordinates for visualisation purposes.
{
    // Obtain all the vector data and then write them as CSV files for visualisation
    std::ofstream file(filename);

    if (!file.is_open())
    {
        cout << " Error in Openning the file " << filename << endl;
        cout << " File : Particle.C :: Function : OutputFile()" << endl;
        cout << " exiting Process()" << endl;
        exit(0);
    }

    file << "x"
         << ","
         << "y"
         << ","
         << "z"
         << ","
         << "deposition"
         << ","
         << "error"
         << ","
         << "escaped"
         << ","
         << "stagnant"
         << ","
         << "cell"
         << ","
         << "prevCell"
         << "\n";
    for (int i = 0; i < position_X.size(); i++)
    {
        int depostionStatus = 0;
        if (isParticleDeposited[i])
            depostionStatus = 1;

        file << position_X[i] << "," << position_Y[i] << "," << position_Z[i] << "," << depostionStatus << "," << isErrorParticle[i] << "," << isEscapedParticle[i] << "," << isStagnantParticle[i] << "," << currentCell[i] << "," << previousCell[i] << "\n";
    }

    file.close();
}

// Things not being considered: validity of the file, number of particles, ordering of the columns
int TParticles::UpdateParticleDetailsFromFile(std::string filename)
{
		std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

		std::string line;
		int counter = -1;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string x_str, y_str, z_str, deposition_str, error_str, escaped_str, stagnant_str, cell_str, prevCell_str;

				if (counter == -1) {
					counter++;
					continue;
				}

        if (std::getline(iss, x_str, ',') &&
            std::getline(iss, y_str, ',') &&
            std::getline(iss, z_str, ',') &&
            std::getline(iss, deposition_str, ',') &&
            std::getline(iss, error_str, ',') &&
            std::getline(iss, escaped_str, ',') &&
            std::getline(iss, stagnant_str, ',') &&
            std::getline(iss, cell_str, ',') &&
            std::getline(iss, prevCell_str)) {
            try {
								position_X[counter] = std::stod(x_str);
								position_Y[counter] = std::stod(y_str);
								position_Z[counter] = std::stod(z_str);
								isParticleDeposited[counter] = std::stoi(deposition_str) == 1 ? true : false;
								isErrorParticle[counter] = std::stoi(error_str) == 1 ? true : false;
								isEscapedParticle[counter] = std::stoi(escaped_str) == 1 ? true : false;
								isStagnantParticle[counter] = std::stoi(stagnant_str) == 1 ? true : false;
								currentCell[counter] = std::stoi(cell_str);
								previousCell[counter] = std::stoi(prevCell_str);
                counter++;
            } catch (const std::exception& e) {
                std::cerr << "Error parsing data: " << e.what() << std::endl;
            }
        } else {
            std::cerr << "Invalid line format: " << line << std::endl;
        }
    }

    file.close();

    std::cout << "Data points read: " << counter << std::endl;

		printUpdatedParticleDetailStats();

		m_ParticlesReleased = counter;
		cout << " Number of particles released : " << m_ParticlesReleased << endl;

		// get substring containing the time from the filename
		std::string time = filename.substr(filename.find_last_of("_") + 1);
		time = time.substr(0, time.find("."));
		// TDatabase::TimeDB->CURRENTTIME = std::stoi(time);

		return std::stoi(time);
}

void TParticles::printUpdatedParticleDetailStats()
{
		int errorParticles = 0;
		int escapedParticles = 0;
		int stagnantParticles = 0;
		int depositedParticles = 0;
		double normX = 0;
		double normY = 0;
		double normZ = 0;

		for (int i = 0; i < position_X.size(); i++) {
			if (isErrorParticle[i])
				errorParticles++;
			if (isEscapedParticle[i])
				escapedParticles++;
			if (isStagnantParticle[i])
				stagnantParticles++;
			if (isParticleDeposited[i])
				depositedParticles++;
			normX += position_X[i] * position_X[i];
			normY += position_Y[i] * position_Y[i];
			normZ += position_Z[i] * position_Z[i];
		}

		cout << " Number of error particles : " << errorParticles << endl;
		cout << " Number of escaped particles : " << escapedParticles << endl;
		cout << " Number of stagnant particles : " << stagnantParticles << endl;
		cout << " Number of deposited particles : " << depositedParticles << endl;
		cout << " Norm of x : " << normX << endl;
		cout << " Norm of y : " << normY << endl;
		cout << " Norm of z : " << normZ << endl;
}

// Checks if a particle is stagnant based on the distance between previous position and current position
bool TParticles::isStagnant(int i) {
		double distance = sqrt(pow(position_X[i] - previousPosition_X[i], 2) + pow(position_Y[i] - previousPosition_Y[i], 2) + pow(position_Z[i] - previousPosition_Z[i], 2));
		if (distance < 0.0001)
				return true;
		else
				return false;
}

// Mark particles as stagnant if they stay in the same area for a long time
void TParticles::detectStagnantParticles() {
		if (m_ParticlesReleased == N_Particles) {
				for (int i = 0; i < N_Particles; i++) {
						if (isParticleDeposited[i] != 1 && isStagnant(i)) {
								isParticleDeposited[i] = true;
								isStagnantParticle[i] = 1;
								m_StagnantParticlesCount++;
						}
						previousPosition_X[i] = position_X[i];
						previousPosition_Y[i] = position_Y[i];
						previousPosition_Z[i] = position_Z[i];
				}
		}
}

// Helper Function
// Function to Interpolate the Velocity of the Particle in next time Step
void TParticles::interpolateNewVelocity(double timeStep, TFEVectFunct3D *VelocityFEVectFunction, TFESpace3D* fespace)
{
    // Constants required for computation
    double densityFluid = 1.1385;
    double densityParticle = 914; // earlier 1266
    double g_x = 0;
    double g_y = 0;
    // double g_z = 1.0 / 51.6414;
    double g_z = 9.81;

    double dynamicViscosityFluid = 0.00001893;
    double lambda = 0.00000007;

    // Particle Diameter = 8 Micrometers
    double particleDiameter = 2.5e-6; // earlier 4e-6

    double mass_particle = (Pi * pow(particleDiameter, 3) * densityParticle) / 6;
    double mass_div_dia = mass_particle / particleDiameter;
    // For the First Term
    double intertialConstant = (3. / 4.) * (densityFluid / densityParticle) * (1 / particleDiameter);

    // For the second term
    double gForceConst_x = g_x;
    double gForceConst_y = g_y;
    double gForceConst_z = g_z;

    int MaxLen;
    int N_Joints;
    const int *TmpLen;
    const int *TmpFV;

    int ErrorParticles = 0;

    int FirstTime = 1;

    // Lambda Function to compute CD/CC
    auto CD_CC = [&](double particleVel, double fluidVel)
    {
        double Re_Particle = densityFluid * particleDiameter * fabs(fluidVel - particleVel) / dynamicViscosityFluid;
        double CD = (24 / Re_Particle) * (1 + 0.15 * pow(Re_Particle, 0.687));
        double CC = 1.0 + ((2 * lambda) / particleDiameter) * (1.257 + 0.4 * exp(-1.0 * ((1.1 * particleDiameter) / (2 * lambda))));
        CC = 1.0;
        return CD/CC;
        // return 1;
    };

    // Here We ensure that the Particles are released in timely manner , in batches of 2000, every 10 time steps
    int numParticlesReleasedPerTimeStep = 500;
    int timeStepCounter = 0;
    int timeStepInterval = 10;   // Release particles every n steps
    
    int actualTimeStep = (int) (TDatabase::TimeDB->CURRENTTIME / TDatabase::TimeDB->TIMESTEPLENGTH);


    // release at first time step and at every 10th time step
    if(actualTimeStep % timeStepInterval == 0  || (m_ParticlesReleased ==0))
    {
        m_ParticlesReleased += numParticlesReleasedPerTimeStep;
        cout << " Addional Particles Released : " << numParticlesReleasedPerTimeStep << " Total Particles Released : " << m_ParticlesReleased <<endl;
    }

    if(m_ParticlesReleased > N_Particles)
    {
        m_ParticlesReleased = N_Particles;
        cout << " All Particles Released : " << m_ParticlesReleased <<endl;
    }
    
    //Get the FEFunction3D for the Velocity
    TFEFunction3D *FEFuncVelocityX = VelocityFEVectFunction->GetComponent(0);
    TFEFunction3D *FEFuncVelocityY = VelocityFEVectFunction->GetComponent(1);
    TFEFunction3D *FEFuncVelocityZ = VelocityFEVectFunction->GetComponent(2);

    

    for (int i = 0; i < m_ParticlesReleased; i++)
    {
        // cout << " =====================================================================================================================  " <<endl;
        if (isParticleDeposited[i] == true)
            continue;
        double values[4];
        int CellNo = currentCell[i];
        // cout << "Thread : "  << omp_get_thread_num() << " , " << CellNo <<endl;
        
        TBaseCell *cell = fespace->GetCollection()->GetCell(CellNo);
        FEFuncVelocityX->FindValueLocal(cell, CellNo, position_X[i], position_Y[i], position_Z[i], values);
        double fluidVelocityX = values[0];
        
        FEFuncVelocityY->FindValueLocal(cell, CellNo, position_X[i], position_Y[i], position_Z[i], values);
        double fluidVelocityY = values[0];
        FEFuncVelocityZ->FindValueLocal(cell, CellNo, position_X[i], position_Y[i], position_Z[i], values);
        double fluidVelocityZ = values[0];

        cout << "----------------- Particle " << i  <<" : -------------------------------------------" <<endl;
        cout << "Cell No   : " << CellNo <<" ,  " <<currentCell[i]  <<endl;
        cout << "POsition X: " <<  position_X[i] <<endl;
        cout << "POsition Y: " <<  position_Y[i] <<endl;
        cout << "POsition Z: " <<  position_Z[i] <<endl;
        cout << "Velocity X: " <<  fluidVelocityX<<endl;
        cout << "Velocity Y: " <<  fluidVelocityY<<endl;
        cout << "Velocity Z: " <<  fluidVelocityZ<<endl;


        // #pragma omp critical
        // cout << i << "," << CellNo << "," << position_X[i] << "," << position_Y[i] << "," << position_Z[i] << "," << fluidVelocityX << "," << fluidVelocityY << "," << fluidVelocityZ << " , " << omp_get_thread_num() << endl;
        

        double cdcc_x = CD_CC(velocityX[i], fluidVelocityX);
        double cdcc_y = CD_CC(velocityY[i], fluidVelocityY);
        double cdcc_z = CD_CC(velocityZ[i], fluidVelocityZ);

        // The RHS will be
        double rhs_x = intertialConstant * cdcc_x * fabs(fluidVelocityX - velocityX[i]) * (fluidVelocityX - velocityX[i]) + gForceConst_x * (densityFluid - densityParticle) / densityParticle;
        double rhs_y = intertialConstant * cdcc_y * fabs(fluidVelocityY - velocityY[i]) * (fluidVelocityY - velocityY[i]) + gForceConst_y * (densityFluid - densityParticle) / densityParticle;
        double rhs_z = intertialConstant * cdcc_z * fabs(fluidVelocityZ - velocityZ[i]) * (fluidVelocityZ - velocityZ[i]) + gForceConst_z * (densityFluid - densityParticle) / densityParticle;

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
        for (int cellId = 0; cellId < N_Cells; cellId++)
        {
            cell = fespace->GetCollection()->GetCell(cellId);

            if (cell->PointInCell(position_X[i], position_Y[i], position_Z[i]))
            {
                insideDomain = true;
                previousCell[i] = currentCell[i];
                currentCell[i] = cellId;
                break;
            }
        }

        cout << "Particle " << i << " is in cell " << currentCell[i] << " and was in cell " << previousCell[i] << endl;

        // If Not, particle is deposited.
        if (!insideDomain)
        {
            // Find the last Cell of the particle
            int cellNo = currentCell[i];
            // int  jointID= Face_id_cellsOnBoundary[i];

            TBaseCell *cell = fespace->GetCollection()->GetCell(cellNo);

            FE3D elementId = fespace->GetFE3D(cellNo, cell);
            TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
            TFEDesc3D *fedesc = element->GetFEDesc3D();
            cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

            int **JointDOFs = fedesc->GetJointDOF();

            TBoundFace *Bdface; // Pointer to Boundary Face in 3D Cell
            TBoundComp *BoundComp;
            int N_Joints = cell->GetN_Joints();

            // Identify the Boundary Facedd
            int jointID = -999999999;
            for (int jj = 0; jj < N_Joints; jj++)
            {
                TJoint *Joint = cell->GetJoint(jj);

                // Pick the cells with Boundary Surfaces and Ignore Inlet Boundary faces.
                if (Joint->GetType() == BoundaryFace)
                {
                    Bdface = (TBoundFace *)Joint;
                    BoundComp = Bdface->GetBoundComp();
                    int bdid = BoundComp->GetID();
                    if (bdid != 0)
                    {
                        jointID = jj;
                        break;
                    }
                }
            }

            if (jointID == -999999999) // Last cell was not a Boundary cell
            {
                cout <<"Uncertain Case "<<endl;
                // Check if any DOF within the Cell is dirichlet and update that as deposition.
                TBaseCell *cell = fespace->GetCollection()->GetCell(currentCell[i]);

                int *globalDOFindex = fespace->GetGlobalNumbers();
                int *BeginIndex = fespace->GetBeginIndex();
                int start = BeginIndex[currentCell[i]];
                int end = BeginIndex[currentCell[i] + 1];
                int N_Active = fespace->GetActiveBound();

                int boundaryDOFfound = 0;
                int foundDOF = 0;
                for (int index = start; index < end; index++)
                {
                    if (globalDOFindex[index] >= N_Active)
                    {
                        boundaryDOFfound = 1;
                        foundDOF = globalDOFindex[index];
                        break;
                    }
                }

                if (boundaryDOFfound == 1)
                {
                    cout <<"Boundary DOF "<<endl;
                    double xx = 0.;
                    double yy = 0.;
                    double zz = 0.;
                    // mark the deposition as the vertex point.
                    fespace->GetDOFPosition(foundDOF, xx, yy, zz);
                    position_X[i] = xx;
                    position_Y[i] = yy;
                    position_Z[i] = zz;
                    isParticleDeposited[i] = true;
                    continue;
                }
                else
                {
                    // Check for Boundary of the neighbours for particles.
                    cout << "Error" << endl;
                    isParticleDeposited[i] = true;
                    m_ErrorParticlesCount++;
                    isErrorParticle[i] = 1;
                    continue;
                }
            }

            TJoint *Joint = cell->GetJoint(jointID);
            double x1, x2, x3, y1, y2, y3, z1, z2, z3;

            cell->GetVertex(TmpFV[jointID * MaxLen + 0])->GetCoords(x1, y1, z1);
            cell->GetVertex(TmpFV[jointID * MaxLen + 1])->GetCoords(x2, y2, z2);
            double t11 = x2 - x1;
            double t12 = y2 - y1;
            double t13 = z2 - z1;
            double len = sqrt(t11 * t11 + t12 * t12 + t13 * t13);
            t11 /= len;
            t12 /= len;
            t13 /= len;

            cell->GetVertex(TmpFV[jointID * MaxLen + (TmpLen[jointID] - 1)])->GetCoords(x2, y2, z2);
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

            VectorStruct firstPoint;
            VectorStruct secondPoint;
            VectorStruct LineVector;
            VectorStruct pointOnSurface;
            VectorStruct normalSurface;
            VectorStruct temp1;
            VectorStruct temp2;

            firstPoint.x = position_X_old[i];
            firstPoint.y = position_Y_old[i];
            firstPoint.z = position_Z_old[i];
            secondPoint.x = position_X[i];
            secondPoint.y = position_Y[i];
            secondPoint.z = position_Z[i];

            normalSurface.x = N1;
            normalSurface.y = N2;
            normalSurface.z = N3;
            pointOnSurface.x = x1;
            pointOnSurface.y = y1;
            pointOnSurface.z = z1;

            //  u = p1 - p0
            LineVector.x = firstPoint.x - secondPoint.x;
            LineVector.y = firstPoint.y - secondPoint.y;
            LineVector.z = firstPoint.z - secondPoint.z;

            // Dot
            double dot = normalSurface.x * firstPoint.x + normalSurface.y * firstPoint.y + normalSurface.z * firstPoint.z;

            if (fabs(dot - 0.0) > 1e-3)
            {
                // w = p0 - pC0
                temp1.x = firstPoint.x - pointOnSurface.x;
                temp1.y = firstPoint.y - pointOnSurface.y;
                temp1.z = firstPoint.z - pointOnSurface.z;

                double fac = -1.0 * (normalSurface.x * temp1.x + normalSurface.y * temp1.y + normalSurface.z * temp1.z);
                fac /= fac;

                // u = u*fac
                LineVector.x = LineVector.x * fac;
                LineVector.y = LineVector.y * fac;
                LineVector.z = LineVector.z * fac;

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

    cout << " Difference in position X : " << Ddot(N_Particles, position_X.data(), position_X_old.data()) << "\n";
    cout << " Difference in position Y : " << Ddot(N_Particles, position_Y.data(), position_Y_old.data()) << "\n";
    cout << " Difference in position Z : " << Ddot(N_Particles, position_Z.data(), position_Z_old.data()) << "\n";

    int depositedCount = 0;
    for (int l = 0; l < isParticleDeposited.size(); l++)
        if (isParticleDeposited[l] == true)
            depositedCount++;
    int NotDeposited = N_Particles - depositedCount;
    cout << "No of particles Deposited : " << depositedCount << endl;
    cout << "percentage of particles Not Deposited : " << (double(N_Particles - depositedCount) / (double)N_Particles) * 100 << " % " << endl;
    cout << "Error particles Accumulaated : " << m_ErrorParticlesCount << endl;
}



// Function to Interpolate the velocity of the particles in next time step 
// PARALLEL Version 
// Function to Interpolate the Velocity of the Particle in next time Step
void TParticles::interpolateNewVelocity_Parallel(double timeStep, TFEVectFunct3D *VelocityFEVectFunction, TFESpace3D* fespace)
{
    // Constants required for computation
    double densityFluid = 1.1385;
    double densityParticle = 914; // earlier 1266
    double g_x = 0;
    double g_y = 0;
    // double g_z = 1.0 / 51.6414;
    double g_z = 9.81;

    double dynamicViscosityFluid = 0.00001893;
    double lambda = 0.00000007;

    // Particle Diameter = 8 Micrometers
    double particleDiameter = 4.3e-6; // earlier 4e-6

    double mass_particle = (Pi * pow(particleDiameter, 3) * densityParticle) / 6;
    double mass_div_dia = mass_particle / particleDiameter;
    // For the First Term
    double intertialConstant = (3. / 4.) * (densityFluid / densityParticle) * (1 / particleDiameter);

    // For the second term
    double gForceConst_x = g_x;
    double gForceConst_y = g_y;
    double gForceConst_z = g_z;

    int MaxLen;
    int N_Joints;
    const int *TmpLen;
    const int *TmpFV;

    int ErrorParticles = 0;

    int FirstTime = 1;

    // Lambda Function to compute CD/CC
    auto CD_CC = [&](double particleVel, double fluidVel)
    {
        double Re_Particle = densityFluid * particleDiameter * fabs(fluidVel - particleVel) / dynamicViscosityFluid;
        double CD = (24 / Re_Particle) * (1 + 0.15 * pow(Re_Particle, 0.687));
        double CC = 1.0 + ((2 * lambda) / particleDiameter) * (1.257 + 0.4 * exp(-1.0 * ((1.1 * particleDiameter) / (2 * lambda))));
        // CC = 1.0;
        return CD/CC;
        // return 1;
    };

    // Here We ensure that the Particles are released in timely manner , in batches of 2000, every 10 time steps
    int numParticlesReleasedPerTimeStep = 1000;
    int timeStepCounter = 0;
    int timeStepInterval = 10;   // Release particles every n steps
    
    int actualTimeStep = (int) (TDatabase::TimeDB->CURRENTTIME / TDatabase::TimeDB->TIMESTEPLENGTH);


    // release at first time step and at every 10th time step
    if(actualTimeStep % timeStepInterval == 0  || (m_ParticlesReleased ==0))
    {
        m_ParticlesReleased += numParticlesReleasedPerTimeStep;
        cout << " Addional Particles Released : " << numParticlesReleasedPerTimeStep << " Total Particles Released : " << m_ParticlesReleased <<endl;
    }

    if(m_ParticlesReleased > N_Particles)
    {
        m_ParticlesReleased = N_Particles;
        cout << " All Particles Released : " << m_ParticlesReleased <<endl;
    }
    
    //Get the FEFunction3D for the Velocity
    TFEFunction3D *FEFuncVelocityX = VelocityFEVectFunction->GetComponent(0);
    TFEFunction3D *FEFuncVelocityY = VelocityFEVectFunction->GetComponent(1);
    TFEFunction3D *FEFuncVelocityZ = VelocityFEVectFunction->GetComponent(2);


    cout << "PART RELWEASED : " << m_ParticlesReleased <<endl;

		int num_threads = (int) ceil(0.9 * omp_get_max_threads());
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < m_ParticlesReleased; i++)
    {
        cout << " =====================================================================================================================  " <<endl;
        if (isParticleDeposited[i] == true)
            continue;
        double values[4];
        int CellNo = currentCell[i];
        // cout << "Thread : "  << omp_get_thread_num() << " , " << CellNo <<endl;
        
        TBaseCell *cell = fespace->GetCollection()->GetCell(CellNo);
        FEFuncVelocityX->FindValueLocal_Parallel(cell, CellNo, position_X[i], position_Y[i], position_Z[i], values);
        double fluidVelocityX = values[0];
        
        FEFuncVelocityY->FindValueLocal_Parallel(cell, CellNo, position_X[i], position_Y[i], position_Z[i], values);
        double fluidVelocityY = values[0];
        FEFuncVelocityZ->FindValueLocal_Parallel(cell, CellNo, position_X[i], position_Y[i], position_Z[i], values);
        double fluidVelocityZ = values[0];

        # pragma omp critical
        {
            cout << "----------------- Particle " << i  <<" : -------------------------------------------" <<endl;
            cout << "Thread : "  << omp_get_thread_num() << " , " << CellNo <<endl;
            cout << "Cell No   : " << CellNo <<endl;
            cout << "POsition X: " <<  position_X[i] <<endl;
            cout << "POsition Y: " <<  position_Y[i] <<endl;
            cout << "POsition Z: " <<  position_Z[i] <<endl;
            cout << "Velocity X: " <<  fluidVelocityX<<endl;
            cout << "Velocity Y: " <<  fluidVelocityY<<endl;
            cout << "Velocity Z: " <<  fluidVelocityZ<<endl;
        }


        #pragma omp critical
        cout << i << "," << CellNo << "," << position_X[i] << "," << position_Y[i] << "," << position_Z[i] << "," << fluidVelocityX << "," << fluidVelocityY << "," << fluidVelocityZ << " , " << omp_get_thread_num() << endl;
        
        exit(0);
        double cdcc_x = CD_CC(velocityX[i], fluidVelocityX);
        double cdcc_y = CD_CC(velocityY[i], fluidVelocityY);
        double cdcc_z = CD_CC(velocityZ[i], fluidVelocityZ);

        // equivalent to setting Re_p as L_inf norm
        cdcc_x = std::min({cdcc_x, cdcc_y, cdcc_z});
        cdcc_y = cdcc_x;
        cdcc_z = cdcc_x;


        // The RHS will be                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        double rhs_x = intertialConstant * cdcc_x * fabs(fluidVelocityX - velocityX[i]) * (fluidVelocityX - velocityX[i]) + gForceConst_x * (densityFluid - densityParticle) / densityParticle;
        double rhs_y = intertialConstant * cdcc_y * fabs(fluidVelocityY - velocityY[i]) * (fluidVelocityY - velocityY[i]) + gForceConst_y * (densityFluid - densityParticle) / densityParticle;
        double rhs_z = intertialConstant * cdcc_z * fabs(fluidVelocityZ - velocityZ[i]) * (fluidVelocityZ - velocityZ[i]) + gForceConst_z * (densityFluid - densityParticle) / densityParticle;

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

   

        // Update the current position of cell.
        // Check if the particle Exists in the current Domain
        int N_Cells = fespace->GetN_Cells();
        bool insideDomain = false;
        for (int cellId = 0; cellId < N_Cells; cellId++)
        {
            TBaseCell* cell = fespace->GetCollection()->GetCell(cellId);
            bool insideCell = cell->PointInCell_Parallel(position_X[i], position_Y[i], position_Z[i]);
            if (insideCell)
            {
                // #pragma omp critical
                // {
                //     cout << "--" <<" For particle : " << i <<" Thread : " << omp_get_thread_num() << " , PrevCell " << previousCell[i] << " , Curr cell " << currentCell[i] << " Cell iD : " << cellId <<endl;
                //     cout <<  "--" <<" For particle : " << i <<" Thread : " << omp_get_thread_num() << position_X[i] << " , " << position_Y[i] << " , " << position_Z[i] <<endl;
                //     double x0,y0,z0;
                //     cell->GetVertex(0)->GetCoords(x0, y0, z0);
                //     cout <<  "--" <<" For particle : " << i <<" Thread : " << omp_get_thread_num() << " x0: " << x0 << " , y0: " << y0 << " , z0: " << z0 <<endl;
                //     double x1,y1,z1;
                //     cell->GetVertex(1)->GetCoords(x1, y1, z1);
                //     cout <<  "--" <<" For particle : " << i <<" Thread : " << omp_get_thread_num() << " x1: " << x1 << " , y1: " << y1 << " , z1: " << z1 <<endl;
                //     double x2,y2,z2;
                //     cell->GetVertex(2)->GetCoords(x2, y2, z2);
                //     cout <<  "--" <<" For particle : " << i <<" Thread : " << omp_get_thread_num() << " x2: " << x2 << " , y2: " << y2 << " , z2: " << z2 <<endl;
                //     double x3,y3,z3;
                //     cell->GetVertex(3)->GetCoords(x3, y3, z3);
                //     cout <<  "--" <<" For particle : " << i <<" Thread : " << omp_get_thread_num() << " x3: " << x3 << " , y3: " << y3 << " , z3: " << z3 <<endl;
                // }
                
                insideDomain = true;
                previousCell[i] = currentCell[i];
                currentCell[i] = cellId;
                
                break;
            }
        }
        
        // #pragma omp critical
        // {
        //     cout << "Particle " << i << " is in cell " << currentCell[i] << " and was in cell " << previousCell[i] << endl;
        // }

        // If Not, particle is deposited.
        if (!insideDomain)
        {
            // Find the last Cell of the particle
            int cellNo = currentCell[i];
            TBaseCell* cell = fespace->GetCollection()->GetCell(cellNo);
            cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
            // int  jointID= Face_id_cellsOnBoundary[i];

            // Boolean Boundary Cell
            bool isBoundaryCell = false;
            int cornerID;
            int jointID;
            bool isBoundaryDOFPresent;
            // if the cellNo is present in map then make isBoundaryCell as true
            #pragma omp critical
            {
                if (m_mapBoundaryFaceIds.find(cellNo) != m_mapBoundaryFaceIds.end())
                {
                    isBoundaryCell = true;
                    // Check the corner id of the cell
                    cornerID = m_cornerTypeOfBoundCells[cellNo];                
                }

                isBoundaryDOFPresent = false;

                // if cell no is present in the BoundaryDOF map, then make isBoundaryDOFPresent as true
                if (m_BoundaryDOFsOnCell.find(cellNo) != m_BoundaryDOFsOnCell.end())
                {
                    isBoundaryDOFPresent = true;
                }
            }

            // Its a Boundary Cell and the cornerID is 21
            // Mark it as escaped. 
            if (isBoundaryCell && cornerID == 21) 
            {
                isParticleDeposited[i] = true;
                m_EscapedParticlesCount++;
                isEscapedParticle[i] = 1;
                continue;
            }
            


            // Last cell was not a boundary Cell
            if (!isBoundaryCell) // Last cell was not a Boundary cell
            {

                if (isBoundaryDOFPresent)  // Boundary DOF is present
                {
		    
                    std::vector<double> boundaryDOF;
		    #pragma omp critical
		    boundaryDOF= m_BoundaryDOFsOnCell[cellNo];

                    position_X[i] = boundaryDOF[0];
                    position_Y[i] = boundaryDOF[1];
                    position_Z[i] = boundaryDOF[2];
                    isParticleDeposited[i] = true;
                    continue;
                }
                else
                {
                    // Check for Boundary of the neighbours for particles.
                    cout << "Error" << endl;
                    isParticleDeposited[i] = true;
                    m_ErrorParticlesCount++;
                    isErrorParticle[i] = 1;
                    continue;
                }
            }

             // Check if the particle is from a corner shared by two faces 0 and 1
            if (cornerID == 20)
            {
                jointID  = 2;
            }
            else
            {
                jointID = cornerID;
            }
            
            TJoint *Joint = cell->GetJoint(jointID);
            double x1, x2, x3, y1, y2, y3, z1, z2, z3;

            cell->GetVertex(TmpFV[jointID * MaxLen + 0])->GetCoords(x1, y1, z1);
            cell->GetVertex(TmpFV[jointID * MaxLen + 1])->GetCoords(x2, y2, z2);
            double t11 = x2 - x1;
            double t12 = y2 - y1;
            double t13 = z2 - z1;
            double len = sqrt(t11 * t11 + t12 * t12 + t13 * t13);
            t11 /= len;
            t12 /= len;
            t13 /= len;

            cell->GetVertex(TmpFV[jointID * MaxLen + (TmpLen[jointID] - 1)])->GetCoords(x2, y2, z2);
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

            VectorStruct firstPoint;
            VectorStruct secondPoint;
            VectorStruct LineVector;
            VectorStruct pointOnSurface;
            VectorStruct normalSurface;
            VectorStruct temp1;
            VectorStruct temp2;

            firstPoint.x = position_X_old[i];
            firstPoint.y = position_Y_old[i];
            firstPoint.z = position_Z_old[i];
            secondPoint.x = position_X[i];
            secondPoint.y = position_Y[i];
            secondPoint.z = position_Z[i];

            normalSurface.x = N1;
            normalSurface.y = N2;
            normalSurface.z = N3;
            pointOnSurface.x = x1;
            pointOnSurface.y = y1;
            pointOnSurface.z = z1;

            //  u = p1 - p0
            LineVector.x = firstPoint.x - secondPoint.x;
            LineVector.y = firstPoint.y - secondPoint.y;
            LineVector.z = firstPoint.z - secondPoint.z;

            // Dot
            double dot = normalSurface.x * firstPoint.x + normalSurface.y * firstPoint.y + normalSurface.z * firstPoint.z;

            if (fabs(dot - 0.0) > 1e-3)
            {
                // w = p0 - pC0
                temp1.x = firstPoint.x - pointOnSurface.x;
                temp1.y = firstPoint.y - pointOnSurface.y;
                temp1.z = firstPoint.z - pointOnSurface.z;

                double fac = -1.0 * (normalSurface.x * temp1.x + normalSurface.y * temp1.y + normalSurface.z * temp1.z);
                fac /= fac;

                // u = u*fac
                LineVector.x = LineVector.x * fac;
                LineVector.y = LineVector.y * fac;
                LineVector.z = LineVector.z * fac;

                temp2.x = LineVector.x + firstPoint.x;
                temp2.y = LineVector.y + firstPoint.y;
                temp2.z = LineVector.z + firstPoint.z;

                // Mark the Particle as deposited
                isParticleDeposited[i] = true;

                

                position_X[i] = temp2.x;
                position_Y[i] = temp2.y;
                position_Z[i] = temp2.z;

                // if jointID is 1, then mark the particle as escaped
                if(jointID == 1)
                {
                    isParticleDeposited[i] = true;
                    m_EscapedParticlesCount++;
                    isEscapedParticle[i] = 1;   // Mark the particle as escaped
                    continue;
                }
            }
            else
            {
                // Mark the Particle as deposited
                isParticleDeposited[i] = true;
                // Ghost particle, Make the vertex as deposition
                position_X[i] = x1;
                position_Y[i] = y1;
                position_Z[i] = z1;
                m_ghostParticlesCount++;

                // if jointID is 1, then mark the particle as escaped
                if(jointID == 1)
                {
                    isParticleDeposited[i] = true;
                    m_EscapedParticlesCount++;
                    isEscapedParticle[i] = 1;   // Mark the particle as escaped
                    continue;
                }
            }
        }
        // check if the particle is in any border cell
    }

    cout << " Difference in position X : " << Ddot(N_Particles, position_X.data(), position_X_old.data()) << "\n";
    cout << " Difference in position Y : " << Ddot(N_Particles, position_Y.data(), position_Y_old.data()) << "\n";
    cout << " Difference in position Z : " << Ddot(N_Particles, position_Z.data(), position_Z_old.data()) << "\n";

    int depositedCount = 0;
    for (int l = 0; l < isParticleDeposited.size(); l++)
        if (isParticleDeposited[l] == true)
            depositedCount++;
    int NotDeposited = N_Particles - depositedCount;

    // cout using right and left allignment and with fixed width
    
    std::cout << std::setw(50) << std::left << "Number of Particles deposited or Escaped" << " : " << std::setw(10) << std::right << depositedCount << std::endl;
    std::cout << std::setw(50) << std::left << "Percentage of Particles Not Deposited" << " : " << std::setw(10) << std::right << (double(N_Particles - depositedCount) / (double)N_Particles) * 100 << " % " << std::endl;
    std::cout << std::setw(50) << std::left << "Error particles Accumulated" << " : " << std::setw(10) << std::right << m_ErrorParticlesCount << std::endl;
    std::cout << std::setw(50) << std::left << "Ghost particles Accumulated" << " : " << std::setw(10) << std::right << m_ghostParticlesCount << std::endl;
    std::cout << std::setw(50) << std::left << "Stagnant particles Accumulated" << " : " << std::setw(10) << std::right << m_StagnantParticlesCount << std::endl;

    // cout << "No of particles Deposited or Escaped: " << depositedCount << endl;
    // cout << "percentage of particles Not Deposited : " << (double(N_Particles - depositedCount) / (double)N_Particles) * 100 << " % " << endl;
    // cout << "Error particles Accumulaated : " << m_ErrorParticlesCount << endl;
}
