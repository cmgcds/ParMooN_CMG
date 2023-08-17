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
#include <cstring>
#include <queue>


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
#include <chrono>

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

// Supporting Functions 
// function to compute the difference between two vectors and their dot product
double diff_dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    // compute the difference between the vectors
    std::vector<double> diff(v1.size());
    for (int i = 0; i < v1.size(); i++) {
        diff[i] = v1[i] - v2[i];
    }

    // compute the dot product of the difference vector
    double dot_product = 0.0;
    for (int i = 0; i < diff.size(); i++) {
        dot_product += diff[i] * diff[i];
    }

    return dot_product;
}

// FUnciton to initialise particle diameters and densities based on mass fractions
void TParticles::GenerateParticleDiameterAndDensity(std::vector<double> &particle_diameters, std::vector<double> &particle_density, \
                                            std::vector<double> mass_fraction, std::vector<double> density_array,
                                             double CMAD, double GSD )
{
    // Calculate the geometric mean and standard deviation in natural logarithm space
    double GM_log = std::log(CMAD);
    double GSD_log = std::log(GSD);
    
    // Calculate the count fraction for each density
    std::vector<double> count_fraction(density_array.size());
    double total_count_fraction = 0.0;

    // calulate the sum of density array vector
    double sum_of_density_array = std::accumulate(density_array.begin(), density_array.end(), 0.0);

    for (int i = 0; i < density_array.size(); i++) {
        count_fraction[i] = mass_fraction[i] * (sum_of_density_array / density_array[i]);
        total_count_fraction += count_fraction[i];
    }

    // normalise the count fraction
    for (int i = 0; i < count_fraction.size(); i++) {
        count_fraction[i] /= total_count_fraction;
    }

    // Generate random numbers following normal distribution
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    // Total particles initialised
    int total_particles_initialised = 0;

    // Generate particles for each density
    for (int i = 0; i < density_array.size(); i++) {
        int num_particles_i = std::round(count_fraction[i] * N_Particles);
        for (int j = 0; j < num_particles_i; j++) {

            if(total_particles_initialised == N_Particles) {
                break;
            }
            // Calculate particle diameter using log-normal distribution formula
            double Z_value = distribution(generator);
            double particle_diameter = std::exp(GM_log + Z_value * GSD_log);

            // Add particle diameter to the vector
            particle_diameters.push_back(particle_diameter);

            particle_density.push_back(density_array[i]);

            total_particles_initialised ++;
        }
    }


    // Fill the remaining particles with the last density
    // Sanity check for round off error while trying to convert the count fraction to number of particles for a given density
    while (total_particles_initialised < N_Particles) {
        particle_diameters.push_back(particle_diameters[total_particles_initialised - 1]);
        particle_density.push_back(particle_density[total_particles_initialised - 1]);
        total_particles_initialised ++;
    }

    std::cout << " PARTICLE DENSITY DISTRIBUTION SUMMARY: \n";
    std::cout << " -------------------------------------- \n";
    // Print the summary of the particle diameters and densities
    for (int i = 0; i < count_fraction[i]; i++) {
        cout <<" Density : " << particle_density[i] << " Number of Particles : " << std::round(count_fraction[i] * N_Particles) << endl;
    }

    std::vector<double> mass_per_density(density_array.size(), 0.0);


    // calculate the mass of particles per density
    for (int i = 0; i < particle_density.size() ;i++)
    {
        for (int j = 0 ; j < density_array.size() ; j++)
        {   
            if (particle_density[i] == density_array[j])
            {
                mass_per_density[j] += (4.0 / 3.0) * M_PI * pow(particle_diameters[i] / 2.0, 3) * particle_density[i];
            }
        }
    }

    // calculate the total mass of particles
    double total_mass = std::accumulate(mass_per_density.begin(), mass_per_density.end(), 0.0);

    std::cout << " ------------------------- \n";
    std::cout << " TOTAL MASS OF PARTICLES : " << total_mass << " kg \n";
    std::cout << " ------------------------- \n";

    // calculate the mass fraction of particles per density
    std::vector<double> mass_fraction_per_density(density_array.size());

    for (int i = 0; i < density_array.size(); i++) {
    mass_fraction_per_density[i] = mass_per_density[i] / total_mass;
    std::cout << std::setw(10) << "Density: " << std::setw(10) << density_array[i]
              << std::setw(15) << "Mass Fraction: " << std::setw(10) << mass_fraction_per_density[i]
              << std::setw(10) << "Number: " << std::setw(10) << std::round(count_fraction[i] * N_Particles) << std::endl;
    }

    // Open a csv file and write the count fraction and mass fraction of particles per density
    std::ofstream csv_file;
    csv_file.open("particle_density_distribution.csv");

    csv_file << "density, given_mass_fraction, count_fraction, count_of_particles, mass_fraction_computed, error\n";

    for (int i = 0; i < density_array.size(); i++) {
        csv_file << density_array[i] << "," << mass_fraction[i] << "," << count_fraction[i] << "," << std::round(count_fraction[i] * N_Particles) << "," << mass_fraction_per_density[i] << "," << mass_fraction[i] - mass_fraction_per_density[i] << "\n";
    }

    // close the csv file
    csv_file.close();

}

/* Constructor for Particle Class */
TParticles::TParticles(int N_Particles_, double circle_x, double circle_y, double radius, TFESpace3D *fespace, std::string resumeFile, int *StartNo)
{
    // Initialise particle Stream
    this->N_Particles = N_Particles_;

    // Resize all the vectors

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
    isGhostParticle.resize(N_Particles, 0);


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

    isParticleDeposited.resize(N_Particles, 0);

    // Create N particles
    Initialiseparticles(N_Particles_, circle_x, circle_y, radius, fespace);

    // INitialise particle and fluid Parameters required for the simulation
    InitialiseParticleParameters(N_Particles_);

		if (resumeFile != "") {
			*StartNo = UpdateParticleDetailsFromFile(resumeFile);
		}

    #ifdef _CUDA

    // time the below function
    auto start_time = std::chrono::high_resolution_clock::now();
    SetupCudaDataStructures(fespace);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    cout << "Time taken for SetupCudaDataStructures : " << elapsed.count() << " s\n";
    #endif


}


// Initialise particle and fluid Parameters required for the simulation
void TParticles::InitialiseParticleParameters(int N_Particles_)
{
    // Density of the fluid
    m_fluid_density = 1.1385; // kg/m^3

    // Dynamic Viscosity of the fluid
    m_fluid_dynamic_viscosity = 0.00001893; // Pa.s

    // mean free path
    m_lambda = 0.00000007; // m

    // resize particle diameter
    m_particle_diameter.resize(N_Particles_, 4.3e-6); // earlier 10e-6   // [IMPORTANT] THE VALUE NEEDS TO BE CHANGED BELOW

     // resize density
    m_particle_density.resize(N_Particles_, 914); // earlier 1266        // [IMPORTANT] THE VALUE NEEDS TO BE CHANGED BELOW

    // Provide the mass fraction of the particles
    // f propylene glycol, glycerol, nicotine and water with their mass fractions being 52.6%, 28.1%, 1.9%, and 17.4%
    std::vector<double> mass_fraction = {.526, 0.281, 0.019, 0.174};

     // Propelene Glycon : https://www.chemeo.com/cid/23-447-0/Propylene-Glycol  Temp : 293K
     // glycerol : 1260 https://en.wikipedia.org/wiki/Glycerol
     // nicotine : 1010 https://en.wikipedia.org/wiki/Nicotine
    std::vector<double> density_array = {1032.56, 1261, 1010, 1000};

    // Set all the particle sizes bsed on the function 
    std::string particle_size_distribution = "polydisperse";  // or "monodisperse"

    // Set the polydisperse type
    // Fused will generate all particles with same density equivalent to the density difference of all particles. 
    // Seperate will generate particles with different densities
    std::string polydisperse_type = "seperate"; // or "fused"

    // check if the string equals to polydisperse
    if(particle_size_distribution == "polydisperse")
    {
        TParticles::GenerateParticleDiameterAndDensity(m_particle_diameter, m_particle_density, mass_fraction, density_array, 0.307 * 1e-6, 1.69);

        if (polydisperse_type == "fused") {
            double density_sum = 0.0;
            for (int i = 0; i < mass_fraction.size(); i++) {
                density_sum += mass_fraction[i] * density_array[i];
            }

            for (int i = 0; i < m_particle_density.size(); i++) {
                m_particle_density[i] = density_sum;
            }
        }
    }
    else
    {
        // HARDCODED_THIVIN For mono disperse particle with single density 
        double density_mono = 914; // DEHS from sininhale paper
        for (int i = 0; i < m_particle_density.size(); i++) {
            m_particle_density[i] = density_mono;
        }

        // HARDCODED_THIVIN particle size
        double particle_size = 8e-6;
        for (int i = 0; i < m_particle_diameter.size(); i++) {
            m_particle_diameter[i] = particle_size;
        }
    }

    
    
    // --- FOr Polydisperse Fused particle ---- //
    // rewrite all particle density to mass_fraction * density
    

    // // HARDCODED_THIVIN For mono disperse particle with single density 
    // // double density_mono = 1261; // glycerol
    // // for (int i = 0; i < m_particle_density.size(); i++) {
    // //     m_particle_density[i] = density_mono;
    // // }

    // // HARDCODED_THIVIN particle size
    // // double particle_size = 8;
    // for (int i = 0; i < m_particle_diameter.size(); i++) {
    //     m_particle_diameter[i] = particle_size;
    // }


    cout << " INFO : Particle Density has been updated to accumulated value of all densities based on the mass fraction \n";

    // Gravity
    m_gravity_x = 0;
    m_gravity_y = 0;
    m_gravity_z = 9.81; // m/s^2

    // Store the time step 
    #ifdef _CUDA
    h_m_time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    #endif
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

    // Define a vector to store the cells on the boundary
    std::vector<int> cells_on_inlet_boundary;
    int inlet_boundary_id = 0;

    // Lets pick up all the cells, which are at the inlet boundary and store them in a vector
    for(int i = 0 ; i < N_Cells ; i++)
    {
       TBaseCell *cell = fespace->GetCollection()->GetCell(i);

        // NOw Obtain the faces of the cell.
        FE3D elementId = fespace->GetFE3D(i, cell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        // cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        TBoundFace *Bdface; // Pointer to Boundary Face in 3D Cell
        TBoundComp *BoundComp;
        int N_Joints = cell->GetN_Joints();
        bool BoundaryJointFound = false;

        for (int jointId = 0; jointId < N_Joints; jointId++)
        {
            TJoint *Joint = cell->GetJoint(jointId);

            // Pick the cells with Boundary Surfaces and Ignore Inlet Boundary faces.
            if (Joint->GetType() == BoundaryFace)
            {
                Bdface = (TBoundFace *)Joint;
                BoundComp = Bdface->GetBoundComp();
                int bdid = BoundComp->GetID();
                if(bdid == inlet_boundary_id)
                {
                    cells_on_inlet_boundary.push_back(i);
                    continue;
                }
            }
        }

        // get the centroid of the cell and check if it lies within 0.01 from the point (0,0,0)
        double centroid[3];
        double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

        // get the cell vertices
        cell->GetVertex(0)->GetCoords(x0, y0, z0);
        cell->GetVertex(1)->GetCoords(x1, y1, z1);
        cell->GetVertex(2)->GetCoords(x2, y2, z2);
        cell->GetVertex(3)->GetCoords(x3, y3, z3);

        // calculate the centroid
        centroid[0] = (x0 + x1 + x2 + x3) / 4.0;
        centroid[1] = (y0 + y1 + y2 + y3) / 4.0;
        centroid[2] = (z0 + z1 + z2 + z3) / 4.0;

        // check if the centroid lies within 0.01 from the point (0,0,0)
        if (std::sqrt(centroid[0] * centroid[0] + centroid[1] * centroid[1] + centroid[2] * centroid[2]) < 0.01)
        {
            // if the cell lies within 0.01 from the point (0,0,0), then set the cell number as the starting cell
            cells_on_inlet_boundary.push_back(i);
            continue;
        }
    }

    // Print the number of cells on the inlet boundary
    cout << "[INFO] Number of cells on the inlet boundary : " << cells_on_inlet_boundary.size() << endl;

    // write the cells on the inlet boundary to a file
    std::ofstream cells_on_inlet_boundary_file;
    cells_on_inlet_boundary_file.open("cells_on_inlet_boundary.txt");
    
    for (int i = 0; i < cells_on_inlet_boundary.size(); i++) {
        cells_on_inlet_boundary_file << cells_on_inlet_boundary[i] << "\n";
    }

    cells_on_inlet_boundary_file.close();

    // For the siminhale, the inlet in in k0

    for (int particleNo = 0; particleNo < N_Particles; particleNo++)
    {
        if(particleNo % 10000 == 0)
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

//         if(particleNo==0){
// position_X[particleNo]=-0.007369244;
// position_Z[particleNo]=-0.0008269974;
// }
// if(particleNo==1){
// position_X[particleNo]=-0.005620816;
// position_Z[particleNo]=0.003577294;
// }
// if(particleNo==2){
// position_X[particleNo]=0.008693858;
// position_Z[particleNo]=0.0003883274;
// }
// if(particleNo==3){
// position_X[particleNo]=-0.009308558;
// position_Z[particleNo]=0.0005940039;
// }
// if(particleNo==4){
// position_X[particleNo]=0.003735454;
// position_Z[particleNo]=0.00860873;
// }
// if(particleNo==5){
// position_X[particleNo]=0.0005385756;
// position_Z[particleNo]=0.003078379;
// }
// if(particleNo==6){
// position_X[particleNo]=0.004023812;
// position_Z[particleNo]=0.005243961;
// }
// if(particleNo==7){
// position_X[particleNo]=0.00512821;
// position_Z[particleNo]=-0.002693227;
// }
// if(particleNo==8){
// position_X[particleNo]=-0.001271772;
// position_Z[particleNo]=-0.0004453647;
// }
// if(particleNo==9){
// position_X[particleNo]=-0.004501863;
// position_Z[particleNo]=-0.006669856;
// }
// if(particleNo==10){
// position_X[particleNo]=9.04579e-05;
// position_Z[particleNo]=-0.003619341;
// }
// if(particleNo==11){
// position_X[particleNo]=-0.0001204663;
// position_Z[particleNo]=-0.008185342;
// }
// if(particleNo==12){
// position_X[particleNo]=-0.008525018;
// position_Z[particleNo]=-0.002317157;
// }
// if(particleNo==13){
// position_X[particleNo]=0.008276349;
// position_Z[particleNo]=-0.0007110835;
// }
// if(particleNo==14){
// position_X[particleNo]=-0.007492692;
// position_Z[particleNo]=0.003769106;
// }
// if(particleNo==15){
// position_X[particleNo]=0.002590868;
// position_Z[particleNo]=0.00450824;
// }
// if(particleNo==16){
// position_X[particleNo]=0.007771444;
// position_Z[particleNo]=-0.003873563;
// }
// if(particleNo==17){
// position_X[particleNo]=0.000265474;
// position_Z[particleNo]=0.006919631;
// }
// if(particleNo==18){
// position_X[particleNo]=0.006830213;
// position_Z[particleNo]=-0.001692108;
// }
// if(particleNo==19){
// position_X[particleNo]=-0.0006416526;
// position_Z[particleNo]=-0.006433446;
// }
// if(particleNo==20){
// position_X[particleNo]=0.001433096;
// position_Z[particleNo]=-0.009338925;
// }
// if(particleNo==21){
// position_X[particleNo]=-3.039762e-05;
// position_Z[particleNo]=0.004965853;
// }
// if(particleNo==22){
// position_X[particleNo]=-0.00574497;
// position_Z[particleNo]=-0.007391455;
// }
// if(particleNo==23){
// position_X[particleNo]=-0.004508237;
// position_Z[particleNo]=-0.001714135;
// }
// if(particleNo==24){
// position_X[particleNo]=0.004196392;
// position_Z[particleNo]=-0.005201784;
// }
// if(particleNo==25){
// position_X[particleNo]=-0.003649209;
// position_Z[particleNo]=0.003041174;
// }
// if(particleNo==26){
// position_X[particleNo]=0.003626924;
// position_Z[particleNo]=-0.002245493;
// }
// if(particleNo==27){
// position_X[particleNo]=-0.001824666;
// position_Z[particleNo]=0.001297974;
// }
// if(particleNo==28){
// position_X[particleNo]=-0.000229709;
// position_Z[particleNo]=0.009221903;
// }
// if(particleNo==29){
// position_X[particleNo]=-0.006004856;
// position_Z[particleNo]=0.002585383;
// }
// if(particleNo==30){
// position_X[particleNo]=0.003025075;
// position_Z[particleNo]=0.00606146;
// }
// if(particleNo==31){
// position_X[particleNo]=-0.0004713639;
// position_Z[particleNo]=-0.005934993;
// }
// if(particleNo==32){
// position_X[particleNo]=-0.001793739;
// position_Z[particleNo]=0.007712967;
// }
// if(particleNo==33){
// position_X[particleNo]=-0.006756029;
// position_Z[particleNo]=-0.002693219;
// }
// if(particleNo==34){
// position_X[particleNo]=-0.007297813;
// position_Z[particleNo]=-0.0008938543;
// }
// if(particleNo==35){
// position_X[particleNo]=-0.0009539966;
// position_Z[particleNo]=0.008633488;
// }
// if(particleNo==36){
// position_X[particleNo]=0.007217197;
// position_Z[particleNo]=0.0001191175;
// }
// if(particleNo==37){
// position_X[particleNo]=0.00635123;
// position_Z[particleNo]=-0.0007551004;
// }
// if(particleNo==38){
// position_X[particleNo]=0.002654774;
// position_Z[particleNo]=0.006493948;
// }
// if(particleNo==39){
// position_X[particleNo]=-0.004213676;
// position_Z[particleNo]=0.0002886933;
// }
// if(particleNo==40){
// position_X[particleNo]=-0.00171943;
// position_Z[particleNo]=0.007531314;
// }
// if(particleNo==41){
// position_X[particleNo]=0.004594954;
// position_Z[particleNo]=0.004312847;
// }
// if(particleNo==42){
// position_X[particleNo]=0.0004997481;
// position_Z[particleNo]=-0.008696123;
// }
// if(particleNo==43){
// position_X[particleNo]=-0.0002211368;
// position_Z[particleNo]=0.003640982;
// }
// if(particleNo==44){
// position_X[particleNo]=-0.001079531;
// position_Z[particleNo]=0.0002931892;
// }
// if(particleNo==45){
// position_X[particleNo]=-0.001205489;
// position_Z[particleNo]=0.006132998;
// }
// if(particleNo==46){
// position_X[particleNo]=-0.005769621;
// position_Z[particleNo]=-0.00692791;
// }
// if(particleNo==47){
// position_X[particleNo]=0.0045467;
// position_Z[particleNo]=-0.001645526;
// }
// if(particleNo==48){
// position_X[particleNo]=0.00361124;
// position_Z[particleNo]=0.006728398;
// }
// if(particleNo==49){
// position_X[particleNo]=0.002591437;
// position_Z[particleNo]=-0.005729064;
// }
// if(particleNo==50){
// position_X[particleNo]=-0.002223535;
// position_Z[particleNo]=0.008950901;
// }
// if(particleNo==51){
// position_X[particleNo]=-0.004615705;
// position_Z[particleNo]=-0.004319299;
// }
// if(particleNo==52){
// position_X[particleNo]=0.005677304;
// position_Z[particleNo]=-0.004356882;
// }
// if(particleNo==53){
// position_X[particleNo]=0.006394522;
// position_Z[particleNo]=-0.002037126;
// }
// if(particleNo==54){
// position_X[particleNo]=-0.006462393;
// position_Z[particleNo]=-0.006845376;
// }
// if(particleNo==55){
// position_X[particleNo]=-0.004856627;
// position_Z[particleNo]=-0.007967252;
// }
// if(particleNo==56){
// position_X[particleNo]=0.002694349;
// position_Z[particleNo]=0.005895396;
// }
// if(particleNo==57){
// position_X[particleNo]=0.005058808;
// position_Z[particleNo]=0.002668598;
// }
// if(particleNo==58){
// position_X[particleNo]=0.001964333;
// position_Z[particleNo]=-0.003624442;
// }
// if(particleNo==59){
// position_X[particleNo]=-0.007651261;
// position_Z[particleNo]=0.0005224655;
// }
// if(particleNo==60){
// position_X[particleNo]=0.001759778;
// position_Z[particleNo]=0.004059788;
// }
// if(particleNo==61){
// position_X[particleNo]=0.001136716;
// position_Z[particleNo]=-0.003679853;
// }
// if(particleNo==62){
// position_X[particleNo]=0.0005709654;
// position_Z[particleNo]=0.001762382;
// }
// if(particleNo==63){
// position_X[particleNo]=-0.001383044;
// position_Z[particleNo]=-0.002595474;
// }
// if(particleNo==64){
// position_X[particleNo]=-0.001062681;
// position_Z[particleNo]=-0.002243376;
// }
// if(particleNo==65){
// position_X[particleNo]=0.003524738;
// position_Z[particleNo]=0.004572167;
// }
// if(particleNo==66){
// position_X[particleNo]=-0.0007913165;
// position_Z[particleNo]=0.00322711;
// }
// if(particleNo==67){
// position_X[particleNo]=0.002112794;
// position_Z[particleNo]=-0.006992122;
// }
// if(particleNo==68){
// position_X[particleNo]=-0.003123649;
// position_Z[particleNo]=0.0004561543;
// }
// if(particleNo==69){
// position_X[particleNo]=0.006879507;
// position_Z[particleNo]=-0.000884968;
// }
// if(particleNo==70){
// position_X[particleNo]=-0.001267231;
// position_Z[particleNo]=0.005034202;
// }
// if(particleNo==71){
// position_X[particleNo]=0.003913583;
// position_Z[particleNo]=0.005077255;
// }
// if(particleNo==72){
// position_X[particleNo]=0.008066024;
// position_Z[particleNo]=0.004933583;
// }
// if(particleNo==73){
// position_X[particleNo]=0.0008804542;
// position_Z[particleNo]=-0.005683482;
// }
// if(particleNo==74){
// position_X[particleNo]=0.002537227;
// position_Z[particleNo]=0.004186397;
// }
// if(particleNo==75){
// position_X[particleNo]=-0.004394222;
// position_Z[particleNo]=-0.002561806;
// }
// if(particleNo==76){
// position_X[particleNo]=-0.003796621;
// position_Z[particleNo]=-0.007924054;
// }
// if(particleNo==77){
// position_X[particleNo]=-0.006794522;
// position_Z[particleNo]=-0.004059558;
// }
// if(particleNo==78){
// position_X[particleNo]=-0.0004439137;
// position_Z[particleNo]=0.0003997697;
// }
// if(particleNo==79){
// position_X[particleNo]=0.00155552;
// position_Z[particleNo]=-0.003376501;
// }
// if(particleNo==80){
// position_X[particleNo]=-0.00941289;
// position_Z[particleNo]=-0.0002069659;
// }
// if(particleNo==81){
// position_X[particleNo]=0.004779175;
// position_Z[particleNo]=0.0007253354;
// }
// if(particleNo==82){
// position_X[particleNo]=0.00228251;
// position_Z[particleNo]=-0.007481344;
// }
// if(particleNo==83){
// position_X[particleNo]=-0.004819716;
// position_Z[particleNo]=-0.001513283;
// }
// if(particleNo==84){
// position_X[particleNo]=-0.0009629754;
// position_Z[particleNo]=0.003310069;
// }
// if(particleNo==85){
// position_X[particleNo]=-0.005808754;
// position_Z[particleNo]=0.0009394967;
// }
// if(particleNo==86){
// position_X[particleNo]=-0.006480515;
// position_Z[particleNo]=-0.001361834;
// }
// if(particleNo==87){
// position_X[particleNo]=0.00374829;
// position_Z[particleNo]=0.006724131;
// }
// if(particleNo==88){
// position_X[particleNo]=0.005388748;
// position_Z[particleNo]=-0.004996895;
// }
// if(particleNo==89){
// position_X[particleNo]=0.005382281;
// position_Z[particleNo]=0.004985032;
// }
// if(particleNo==90){
// position_X[particleNo]=-0.005566654;
// position_Z[particleNo]=0.003052438;
// }
// if(particleNo==91){
// position_X[particleNo]=-0.002570685;
// position_Z[particleNo]=-0.00532617;
// }
// if(particleNo==92){
// position_X[particleNo]=0.0003834973;
// position_Z[particleNo]=-0.008692875;
// }
// if(particleNo==93){
// position_X[particleNo]=0.004865479;
// position_Z[particleNo]=-0.006572413;
// }
// if(particleNo==94){
// position_X[particleNo]=-0.007896303;
// position_Z[particleNo]=-0.002628757;
// }
// if(particleNo==95){
// position_X[particleNo]=-0.003455449;
// position_Z[particleNo]=0.001535517;
// }
// if(particleNo==96){
// position_X[particleNo]=0.005164915;
// position_Z[particleNo]=0.003995822;
// }
// if(particleNo==97){
// position_X[particleNo]=0.008704832;
// position_Z[particleNo]=-0.002758327;
// }
// if(particleNo==98){
// position_X[particleNo]=-0.0008985522;
// position_Z[particleNo]=-0.00519119;
// }
// if(particleNo==99){
// position_X[particleNo]=-0.001568427;
// position_Z[particleNo]=-0.004204671;
// }

        position_X[particleNo] = y;
        position_Y[particleNo] = 0.001;
        position_Z[particleNo] = z;


        //    cout << "Size : " << cells_on_inlet_boundary.size() << endl;

        // Identify, which cell the particle Belongs to.
        int N_Cells = fespace->GetCollection()->GetN_Cells();
        int num_threads = (int) ceil(0.9 * omp_get_max_threads());
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic,1)  shared(cells_on_inlet_boundary, fespace, particleNo,currentCell) 
        for (int index = 0; index < cells_on_inlet_boundary.size(); index++)
        {
            int cellId = cells_on_inlet_boundary[index];  // Added to not to change any further variables within the code
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

    std::ofstream file("ParticlePositions.txt");
    for(int i = 0; i < N_Particles; i++)
    {
        file << position_X[i] << "," << position_Y[i] << "," << position_Z[i] << "," << currentCell[i] << endl;
    }
    

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
        std::vector<int> joint_ids;
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
                joint_ids.push_back(jointId);
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
                // Assign the joint id , based on the index of the 2 in the tmp_face_ids
                // This is required to identify the deposition of the particle.
                #pragma omp critical
                {
                    auto it = std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2);
                    int index = std::distance(tmp_face_ids.begin(), it);
                    m_jointidOfBoundCells[cellNo] = joint_ids[index];
                }

                N_corner_20++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 1) != tmp_face_ids.end() && std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2) != tmp_face_ids.end())
            {
                cornerType = 21;
                // Assign the joint id , based on the index of the 1 in the tmp_face_ids
                // We will consider this particle as escaping particle via bdid 1

                #pragma omp critical
                {
                    auto it = std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 1);
                    int index = std::distance(tmp_face_ids.begin(), it);
                    m_jointidOfBoundCells[cellNo] = joint_ids[index];
                }
                N_corner_21++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2) != tmp_face_ids.end())
            {
                cornerType = 2;
                #pragma omp critical
                {
                    auto it = std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 2);
                    int index = std::distance(tmp_face_ids.begin(), it);
                    m_jointidOfBoundCells[cellNo] = joint_ids[index];
                }
                N_corner_2++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 0) != tmp_face_ids.end())
            {
                cornerType = 0;
                #pragma omp critical
                {
                    auto it = std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 0);
                    int index = std::distance(tmp_face_ids.begin(), it);
                    m_jointidOfBoundCells[cellNo] = joint_ids[index];
                }
                N_corner_0++;
            }
            else if(std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 1) != tmp_face_ids.end())
            {
                cornerType = 1;
                #pragma omp critical
                {
                    auto it = std::find(tmp_face_ids.begin(), tmp_face_ids.end(), 1);
                    int index = std::distance(tmp_face_ids.begin(), it);
                    m_jointidOfBoundCells[cellNo] = joint_ids[index];
                }
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
            int shift_in_index_per_cell = index - start;
            if (globalDOFindex[index] >= N_Active)
            {
                double xx = 0.;
                double yy = 0.;
                double zz = 0.;
                double *xi, *eta, *zeta;
                int N_Points ;
                double absdetjk[1];
                // mark the deposition as the vertex point.
                // fespace->GetDOFPosition(globalDOFindex[index], xx, yy, zz);
                

                RefTrans3D* RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();
                RefTrans3D RefTrans = RefTransArray[elementId];
                BF3DRefElements RefElement = TFEDatabase3D::GetRefElementFromFE3D(elementId);

                TNodalFunctional3D *nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(elementId);
                nf->GetPointsForAll(N_Points, xi, eta, zeta);

                switch(RefTrans)
                {
                    case TetraAffin:
                    TTetraAffin* ta_rt = new TTetraAffin();
                    ((TTetraAffin *)ta_rt)->SetCell(cell);
                    ((TTetraAffin *)ta_rt)->GetOrigFromRef(1, xi+shift_in_index_per_cell, eta+shift_in_index_per_cell, zeta+shift_in_index_per_cell,
                                                            &xx, &yy, &zz, absdetjk);
                    delete ta_rt;
                    break;

                    // default case - Error needs to be thrown
                    default:
                    cout << " Error in RefTrans3D " << endl;
                    cout << " File : Particle.C :: Function : InitialiseParticles()" << endl;
                    cout << " exiting Process()" << endl;
                    exit(0);
                }
                co_ordinates.emplace_back(xx);
                co_ordinates.emplace_back(yy);
                co_ordinates.emplace_back(zz);
                
                // #pragma omp critical
                // {
                //     printf("%d, %d, %d, %f, %f, %f \n",omp_get_thread_num(), cellNo, globalDOFindex[index], xx, yy, zz);
                // }
                
                // Assign the DOF Co-ordinates to the map with cellNo as key
                #pragma omp critical
                m_BoundaryDOFsOnCell[cellNo] = co_ordinates;

                 break;
            }
        }
    }

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

    // Read the adjacency values of the current mesh
    ReadAdjacencyValues("adjacency_matrix_for_siminhale_reflevel_1.txt",true,row_pointer,col_index,fespace);
}

/* Function to read the adjacency values of the current mesh*/
        // If the Parameter isAdjacencyFile is true, then the function reads the adjacency values from the file
        // else the function reads the adjacency values from the mesh file and generates the adjacency values
void TParticles::ReadAdjacencyValues(std::string filename, bool isAdjacencyFile, std::vector<int> &row_ptr, std::vector<int> &col_ind, TFESpace3D* fespace)
{   
    // If the Given File is an adjacency file, then read the values from the file
    if(isAdjacencyFile)
    {
        // read the adjacency matrix file
        // First line contains the Size of rowPtr and colInd
        std::ifstream file(filename);

        // Check if object is valid
        if(!file)
        {
            std::cerr << "Cannot open the File : "<<filename<<std::endl;
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            exit(0);
        }

        std::string str;
        int line_count = 0;

        int num_cells = fespace->GetCollection()->GetN_Cells();

        // Read the first line
        while (std::getline(file, str) ) {
            line_count++;
            // std::cout  << str << '\n';
            // Check if the line contains the word "Tetrahedra"
            if (line_count == 1) {
                // convert the string to integer
                std::istringstream iss(str);
                
                int row_ptr_size;
                int col_ind_size;
                iss >> row_ptr_size >> col_ind_size;  // Add the number of cells to the variable num_cells

                std::cout << "row_ptr size: " << row_ptr_size << '\n';
                std::cout << "col_ind size: " << col_ind_size << '\n';
                

                // Resize the row_ptr vector
                row_ptr.resize(row_ptr_size);

                // Resize the col_ind vector
                col_ind.resize(col_ind_size);
            }
            // Read the row pointer Array
            if(line_count == 2)
            {
                std::istringstream iss(str);
                for(int i = 0; i < row_ptr.size(); ++i)
                {
                    iss >> row_ptr[i];
                }
            }

            // Read the col_ind Array
            if(line_count == 3)
            {
                std::istringstream iss(str);
                for(int i = 0; i < col_ind.size(); ++i)
                {
                    iss >> col_ind[i];
                }
            }
        }

        // Sanity check
        if (num_cells != row_ptr.size() - 1)
        {
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            std::cout << "Number of cells and row_ptr size mismatch" << std::endl;
            exit(0);
        }

        // Sanity check, Check the last element of row_ptr with size of col_ind
        if (row_ptr[row_ptr.size() - 1] != col_ind.size())
        {
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            std::cout << "Last element of row_ptr " << row_ptr[row_ptr.size() - 1] << "  and col_ind " << col_ind.size()   <<" size mismatch" << std::endl;
            exit(0);
        }
    }

    else
    {
        std::string output_filename = "adjacency_matrix_for_siminhale_reflevel_1.txt";

        //read mesh file
        std::ifstream file(filename);

        // Check if object is valid
        if(!file)
        {
            std::cerr << "Cannot open the File : "<<filename<<std::endl;
            std::cout << "Error on file: " << __FILE__ << " on line: " << __LINE__ << std::endl;
            exit(0);
        }

        std::string str;
        int line_count = 0;
        int num_cells = 0;
        // read untill you reach the word "Tetrahedra"
        // THis will translate the pointer untill the word "Tetrahedra" is reached
        while (std::getline(file, str) ) {
            ++line_count;
            // Check if the line contains the word "Tetrahedra"
            if (str.find("Tetrahedra") != std::string::npos) {
                ++line_count; 
                // read the next line
                std::getline(file, str);
                // convert the string to integer
                std::istringstream iss(str);
                iss >> num_cells;  // Add the number of cells to the variable num_cells
                break;
            }
            
        }

        // std::cout << "Line count: " << line_count << '\n';
        std::cout << "[Adjacency: INFO] Number of cells: " << num_cells << '\n';

        // Read the next num_cells lines and store the vertices of each cell
        std::vector<std::vector<int>> vertices_in_cells(num_cells, std::vector<int>(4));

        for (int i = 0; i < num_cells; ++i) {
            std::getline(file, str);
            std::istringstream iss(str);
            iss >> vertices_in_cells[i][0] >> vertices_in_cells[i][1] >> vertices_in_cells[i][2] >> vertices_in_cells[i][3];
        }

        // Create an adjacency matrix in a CSR format, where each row represents a cell and each column represents a cell that shares a face with the cell in the row

        // resize the row_ptr vector
        row_ptr.resize(num_cells+1);

        std::vector<std::vector<int>> col_ind_for_each_cell(num_cells);

        //initialize row_ptr
        row_ptr[0] = 0;

        // get the total number of omp threads
        int num_threads = omp_get_max_threads();

        // get 90% of the number of threads
        int num_threads_for_calc = (int) (0.9 * num_threads);

        // Loop over all cells
        #pragma omp parallel for num_threads(num_threads_for_calc) schedule(dynamic, 1) \
        shared(row_ptr, col_ind_for_each_cell, vertices_in_cells, num_cells) 
        for (int i = 0; i < num_cells; i++) {
            std::vector<int> col_ind_local;

            if(i%10000 == 0)
                std::cout << "[Adjacency: INFO] Reading Adjacency values from cell: " << i << '\n';

        // Loop over the vertices of the cell
        for (int j = 0; j < 4; j++) {
            // Loop over all cells again
            for (int k = 0; k < num_cells; k++) {
                // Loop over the vertices of the cell
                for (int l = 0; l < 4; l++) {
                    // Check if the cell k shares a face with cell i
                    if (vertices_in_cells[i][j] == vertices_in_cells[k][l] && i != k) {
                        // Check if the cell k is already in the adjacency matrix
                        if (std::find(col_ind_local.begin(), col_ind_local.end(), k) == col_ind_local.end()) {
                            // If not, add it to the adjacency matrix
                            col_ind_local.push_back(k);
                        }
                    }
                }
            }
        }
            row_ptr[i+1] = col_ind_local.size();  // Append the size of the col_ind_local vector to the row_ptr vector
            // Sort the col_ind_local vector
            // std::sort(col_ind_local.begin(), col_ind_local.end());
            // Add the col_ind_local vector to the col_ind vector
            col_ind_for_each_cell[i] = col_ind_local;
        }

        // Loop over the row_ptr vector and add the previous element to the current element to perform a prefix sum
        for (int i = 1; i < row_ptr.size(); ++i) {
            row_ptr[i] += row_ptr[i-1];
        }

        // Loop over the col_ind_for_each_cell vector and add the elements to the col_ind vector
        for (int i = 0; i < col_ind_for_each_cell.size(); ++i) {
            col_ind.insert(col_ind.end(), col_ind_for_each_cell[i].begin(), col_ind_for_each_cell[i].end());
        }
        
        // lets write this to a file 
        std::ofstream outfile;
        outfile.open(output_filename);

        // write the size of the row_ptr and col_ind vectors to the file
        outfile << row_ptr.size() << " " << col_ind.size() << '\n';

        // write the row_ptr vector to the file
        for (int i = 0; i < row_ptr.size(); ++i) {
            outfile << row_ptr[i] << " ";
        }
        outfile << '\n';
        //Write the col_ind vector to the file
        for (int i = 0; i < col_ind.size(); ++i) {
            outfile << col_ind[i] << " ";
        }
    }
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
			   << "vel_x"
         << ","
			   << "vel_y"
         << ","
			   << "vel_z"
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
         << ","
         << "diameter"
         << ","
         << "density"
         << "\n";
    for (int i = 0; i < position_X.size(); i++)
    {
        int depostionStatus = 0;
        if (isParticleDeposited[i])
            depostionStatus = 1;

        file << position_X[i] << ","
						 << position_Y[i] << ","
						 << position_Z[i] << ","
					   << velocityX[i] << ","
					   << velocityY[i] << ","
					   << velocityZ[i] << ","
						 << depostionStatus << ","
						 << isErrorParticle[i] << ","
						 << isEscapedParticle[i] << ","
						 << isStagnantParticle[i] << ","
						 << currentCell[i] << ","
						 << previousCell[i] << ","
						 << m_particle_diameter[i] << ","
						 << m_particle_density[i] << "\n";
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
        std::string x_str, y_str, z_str, vel_x_str, vel_y_str, vel_z_str, deposition_str, error_str, escaped_str, stagnant_str, cell_str, prevCell_str, diameter_str, density_str;

				if (counter == -1) {
					counter++;
					continue;
				}

        if (std::getline(iss, x_str, ',') &&
            std::getline(iss, y_str, ',') &&
            std::getline(iss, z_str, ',') &&
            std::getline(iss, vel_x_str, ',') &&
            std::getline(iss, vel_y_str, ',') &&
            std::getline(iss, vel_z_str, ',') &&
            std::getline(iss, deposition_str, ',') &&
            std::getline(iss, error_str, ',') &&
            std::getline(iss, escaped_str, ',') &&
            std::getline(iss, stagnant_str, ',') &&
            std::getline(iss, cell_str, ',') &&
            std::getline(iss, prevCell_str, ',') &&
            std::getline(iss, diameter_str, ',') &&
            std::getline(iss, density_str)) {
            try {
								position_X[counter] = std::stod(x_str);
								position_Y[counter] = std::stod(y_str);
								position_Z[counter] = std::stod(z_str);
								velocityX[counter] = std::stod(vel_x_str);
								velocityY[counter] = std::stod(vel_y_str);
								velocityZ[counter] = std::stod(vel_z_str);
								isParticleDeposited[counter] = std::stoi(deposition_str);
								isErrorParticle[counter] = std::stoi(error_str) == 1;
								isEscapedParticle[counter] = std::stoi(escaped_str) == 1;
								isStagnantParticle[counter] = std::stoi(stagnant_str) == 1;
								currentCell[counter] = std::stoi(cell_str);
								previousCell[counter] = std::stoi(prevCell_str);
								m_particle_diameter[counter] = std::stod(diameter_str);
								m_particle_density[counter] = std::stod(density_str);
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
		TDatabase::TimeDB->CURRENTTIME = (std::stoi(time) - 1) * TDatabase::TimeDB->TIMESTEPLENGTH;
		return std::stoi(time) + 1;
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
		double distance = sqrt(pow(position_X[i] - previousPosition_X[i], 2) +
													 pow(position_Y[i] - previousPosition_Y[i], 2) +
													 pow(position_Z[i] - previousPosition_Z[i], 2));
		return (distance < 0.0001) ? true : false;
}

// Mark particles as stagnant if they stay in the same area for a long time
void TParticles::detectStagnantParticles() {
		if (m_ParticlesReleased == N_Particles) {
				#ifdef _CUDA
						DetectStagnantParticlesHostWrapper(m_ParticlesReleased);
				#else
						int num_threads = (int) ceil(0.9 * omp_get_max_threads());
						#pragma omp parallel for num_threads(num_threads)
						for (int i = 0; i < N_Particles; i++) {
								if (isParticleDeposited[i] != 1 && isStagnant(i)) {
										isParticleDeposited[i] = 1;
										isStagnantParticle[i] = 1;
								}
								previousPosition_X[i] = position_X[i];
								previousPosition_Y[i] = position_Y[i];
								previousPosition_Z[i] = position_Z[i];
						}
				#endif
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

    // Function send all the velocities 
    

    for (int i = 0; i < m_ParticlesReleased; i++)
    {
        // cout << " =====================================================================================================================  " <<endl;
        if (isParticleDeposited[i] == 1)
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

        cout << "-- intertialConstant : " << intertialConstant   << "  CDCC : " << (CD_CC(velocityX[i],fluidVelocityX)) << " fluidVelocityX: " << fluidVelocityX << "    velocityX : " <<velocityX[i] << "\n";
        cout << "-- intertialConstant     : " << intertialConstant <<"  CDCC : " << (CD_CC(velocityY[i],fluidVelocityY)) <<"     fluidVelocityY: " << fluidVelocityY << "    velocityY : " <<velocityY[i] << "\n";
        cout << "-- intertialConstant : " << intertialConstant << "  CDCC : " << (CD_CC(velocityZ[i],fluidVelocityZ)) <<" fluidVelocityZ: " << fluidVelocityZ << "    velocityZ : " <<velocityZ[i] << "\n";

        // Now, Compute the Updated Velocity of particle Using Forward Euler

        double velocityParticle_X_New = rhs_x * (timeStep) + velocityX[i];
        double velocityParticle_Y_New = rhs_y * (timeStep) + velocityY[i];
        double velocityParticle_Z_New = rhs_z * (timeStep) + velocityZ[i];

        cout << "-- vel new x: " << velocityParticle_X_New << " rhs: " << rhs_x << " t: " << (timeStep) << "\n";
        cout << "-- vel new y: " << velocityParticle_Y_New << " rhs: " << rhs_y << " t : " << (timeStep) << "\n";
        cout << "-- vel new z: " << velocityParticle_Z_New << " rhs: " << rhs_z << " t : " << (timeStep) << "\n";

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
                    isParticleDeposited[i] = 1;
                    continue;
                }
                else
                {
                    // Check for Boundary of the neighbours for particles.
                    cout << "Error" << endl;
                    isParticleDeposited[i] = 1;
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
                isParticleDeposited[i] = 1;

                position_X[i] = temp2.x;
                position_Y[i] = temp2.y;
                position_Z[i] = temp2.z;
            }
            else
            {
                // Mark the Particle as deposited
                isParticleDeposited[i] = 1;
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
        if (isParticleDeposited[l] == 1)
            depositedCount++;
    int NotDeposited = N_Particles - depositedCount;
    cout << "No of particles Deposited : " << depositedCount << endl;
    cout << "percentage of particles Not Deposited : " << (double(N_Particles - depositedCount) / (double)N_Particles) * 100 << " % " << endl;
    cout << "Error particles Accumulaated : " << m_ErrorParticlesCount << endl;
}



// Function to Interpolate the velocity of the particles in next time step 
// PARALLEL Version 
// Function to Interpolate the Velocity of the Particle in next time Step
#ifdef _CUDA
void TParticles::interpolateNewVelocity_Parallel(double timeStep, TFEVectFunct3D *VelocityFEVectFunction, TFESpace3D* fespace)
{

    // Here We ensure that the Particles are released in timely manner , in batches of 2000, every 10 time steps
    int numParticlesReleasedPerTimeStep = 5000;
    int timeStepCounter = 0;
    int timeStepInterval = 100;   // Release particles every n steps
    
    int actualTimeStep = (int) (TDatabase::TimeDB->CURRENTTIME / TDatabase::TimeDB->TIMESTEPLENGTH);


    // release at first time step and at every nth time step
    if( (actualTimeStep % timeStepInterval == 0  || (m_ParticlesReleased ==0) ) &&  (m_ParticlesReleased < N_Particles))
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

    if((actualTimeStep+1) % 100 ==0)
    {
        cout << "PARTICLES CURRENTLY IN DOMAIN : " << m_ParticlesReleased <<endl;
    }
    

    // Call the host wrapper function to call the device kernel
    // This function will be called for every time step
    int n_dof = fespace->GetN_DegreesOfFreedom();
    int n_cells = fespace->GetN_Cells();
    double *fe_function_velocity_x = FEFuncVelocityX->GetValues();
    double *fe_function_velocity_y = FEFuncVelocityY->GetValues();
    double *fe_function_velocity_z = FEFuncVelocityZ->GetValues();


    SetupVelocityValues(fe_function_velocity_x,fe_function_velocity_y,fe_function_velocity_z,m_ParticlesReleased,n_dof);
    
    // Call the kernel function wrapper to perform the interpolation
    // This function will be called for every time step
    // std::cout << "Calling the Interpolate Velocity Host Wrapper" << std::endl;
    InterpolateVelocityHostWrapper(timeStep, m_ParticlesReleased,n_dof,n_cells); 

    

    int depositedCount = std::count(isParticleDeposited.begin(), isParticleDeposited.end(), 1);
    m_ErrorParticlesCount = std::count(isErrorParticle.begin(), isErrorParticle.end(), 1);
    m_ghostParticlesCount = std::count(isGhostParticle.begin(), isGhostParticle.end(), 1);
    m_StagnantParticlesCount = std::count(isStagnantParticle.begin(), isStagnantParticle.end(), 1);

    // Enable printing for the Output console once every 100 timesteps
    if((actualTimeStep+1) % 100 ==0)
    {
        double x_difference_position =  diff_dot_product(position_X,position_X_old);
        double y_difference_position =  diff_dot_product(position_Y,position_Y_old);
        double z_difference_position =  diff_dot_product(position_Z,position_Z_old);
        cout << " Difference in position X : " << x_difference_position<< "\n";
        cout << " Difference in position Y : " << y_difference_position<< "\n";
        cout << " Difference in position Z : " << z_difference_position<< "\n";

        std::cout << std::setw(50) << std::left << "Number of Particles deposited or Escaped" << " : " << std::setw(10) << std::right << depositedCount << std::endl;
        std::cout << std::setw(50) << std::left << "Percentage of Particles Not Deposited" << " : " << std::setw(10) << std::right << (double(N_Particles - depositedCount) / (double)N_Particles) * 100 << " % " << std::endl;
        std::cout << std::setw(50) << std::left << "Error particles Accumulated" << " : " << std::setw(10) << std::right << m_ErrorParticlesCount << std::endl;
        std::cout << std::setw(50) << std::left << "Ghost particles Accumulated" << " : " << std::setw(10) << std::right << m_ghostParticlesCount << std::endl;
        std::cout << std::setw(50) << std::left << "Stagnant particles Accumulated" << " : " << std::setw(10) << std::right << m_StagnantParticlesCount << std::endl;
    }
    


    // cout << "No of particles Deposited or Escaped: " << depositedCount << endl;
    // cout << "percentage of particles Not Deposited : " << (double(N_Particles - depositedCount) / (double)N_Particles) * 100 << " % " << endl;
    // cout << "Error particles Accumulaated : " << m_ErrorParticlesCount << endl;
}
#else
void TParticles::interpolateNewVelocity_Parallel(double timeStep, TFEVectFunct3D *VelocityFEVectFunction, TFESpace3D* fespace)
{

    // time the function
    auto start = std::chrono::high_resolution_clock::now();

    // Constants required for computation
    double densityFluid = 1.1385;
    double g_x = 0;
    double g_y = 0;
    double g_z = 9.81;
    double dynamicViscosityFluid = 0.00001893;
    double lambda = 0.00000007;
    // Here We ensure that the Particles are released in timely manner , in batches of 2000, every 10 time steps
    int numParticlesReleasedPerTimeStep = 50000;
    int timeStepCounter = 0;
    int timeStepInterval = 10;   // Release particles every n steps

    // Search Depth
    int searchDepth = 3;



		auto inertialConstant = [&](double densityParticle, double particleDiameter)
		{ return (3. / 4.) * (densityFluid / densityParticle) * (1 / particleDiameter); };

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
    auto CD_CC = [&](double particleVel, double fluidVel, double particleDiameter)
    {
        double Re_Particle = densityFluid * particleDiameter * fabs(fluidVel - particleVel) / dynamicViscosityFluid;
        double CD = (24 / Re_Particle) * (1 + 0.15 * pow(Re_Particle, 0.687));
        double CC = 1.0 + ((2 * lambda) / particleDiameter) * (1.257 + 0.4 * exp(-1.0 * ((1.1 * particleDiameter) / (2 * lambda))));
        // CC = 1.0;
        return CD/CC;
        // return 1;
    };


    
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


    cout << "PART RELEASED : " << m_ParticlesReleased <<endl;

    
    // Open a file to write the data based on the actual time step
    // std::string filename = "Debug_ParticleData_" + std::to_string(actualTimeStep) + ".csv";
    // std::ofstream myfile;
    // myfile.open(filename);

    int num_threads = (int) ceil(0.9 * omp_get_max_threads());
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < m_ParticlesReleased; i++)
    {
        // cout << " ======================" << i << "=============================== \n" <<endl;
        if (isParticleDeposited[i] == 1)
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


        // # pragma omp critical
        // {
        //     cout << "\n----------------- Particle " << i  <<" : -------------------------------------------\n" <<endl;
        //     cout << "Thread : "  << omp_get_thread_num() << " , " << CellNo <<endl;
        //     cout << "Cell No   : " << CellNo <<endl;
        //     cout << "POsition X: " <<  position_X[i] <<endl;
        //     cout << "POsition Y: " <<  position_Y[i] <<endl;
        //     cout << "POsition Z: " <<  position_Z[i] <<endl;
        //     cout << "Velocity X: " <<  fluidVelocityX<<endl;
        //     cout << "Velocity Y: " <<  fluidVelocityY<<endl;
        //     cout << "Velocity Z: " <<  fluidVelocityZ<<endl;
        // }
        
        double cdcc_x = CD_CC(velocityX[i], fluidVelocityX, m_particle_diameter[i]);
        double cdcc_y = CD_CC(velocityY[i], fluidVelocityY, m_particle_diameter[i]);
        double cdcc_z = CD_CC(velocityZ[i], fluidVelocityZ, m_particle_diameter[i]);

        // // equivalent to setting Re_p as L_inf norm
        // cdcc_x = std::min({cdcc_x, cdcc_y, cdcc_z});
        // cdcc_y = cdcc_x;
        // cdcc_z = cdcc_x;

        // The RHS will be
        double rhs_x = inertialConstant(m_particle_density[i], m_particle_diameter[i]) * cdcc_x
					* fabs(fluidVelocityX - velocityX[i]) * (fluidVelocityX - velocityX[i])
					+ gForceConst_x * (densityFluid - m_particle_density[i]) / m_particle_density[i];
        double rhs_y = inertialConstant(m_particle_density[i], m_particle_diameter[i]) * cdcc_y
					* fabs(fluidVelocityY - velocityY[i]) * (fluidVelocityY - velocityY[i])
					+ gForceConst_y * (densityFluid - m_particle_density[i]) / m_particle_density[i];
        double rhs_z = inertialConstant(m_particle_density[i], m_particle_diameter[i]) * cdcc_z
					* fabs(fluidVelocityZ - velocityZ[i]) * (fluidVelocityZ - velocityZ[i])
					+ gForceConst_z * (densityFluid - m_particle_density[i]) / m_particle_density[i];
        

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

        //   #pragma omp critical
        //   {
        //     myfile << i << "," << CellNo << "," << position_X[i] << "," << position_Y[i] << "," << position_Z[i] << "," << fluidVelocityX << "," << fluidVelocityY << "," << fluidVelocityZ << ", " << cdcc_x << ", " << cdcc_y << ", " << cdcc_z << ", " << rhs_x << ", " << rhs_y << ", " << rhs_z << ", " << velocityX[i] << ", " << velocityY[i] << ", " << velocityZ[i] << ", " << position_X[i] << ", " << position_Y[i] << ", " << position_Z[i] << endl;
        //   }

        // Update the current position of cell.
        // Check if the particle Exists in the current Domain
        int N_Cells = fespace->GetN_Cells();
        bool insideDomain = false;


        bool insideCurrentCell = cell->PointInCell_Parallel(position_X[i], position_Y[i], position_Z[i]);
        if(insideCurrentCell)
        {
            insideDomain = true;
            previousCell[i] = currentCell[i];
            currentCell[i] = CellNo;     // Assign Same cell no
            continue;
        }
        else
        {
            // generate a current queue of cells to search
            std::queue<std::pair<int, int>> cell_queue;
            cell_queue.push({CellNo, 0});   // Push the queue with current cell no and search depth 0
            while(!cell_queue.empty() && !insideDomain)
            {
                int current_cell = cell_queue.front().first;
                int current_depth = cell_queue.front().second;
                cell_queue.pop();

                // Get the current levels neighbours of the cell 
                int start_index = row_pointer[current_cell];
                int end_index = row_pointer[current_cell+1];

                // For cells in the current level
                for(int index = start_index; index < end_index; index++)
                {
                    int neighbourCellNo = col_index[index];
                    TBaseCell* neighbourCell = fespace->GetCollection()->GetCell(neighbourCellNo);
                    bool insideNeighbourCell = neighbourCell->PointInCell_Parallel(position_X[i], position_Y[i], position_Z[i]);
                    if(insideNeighbourCell)
                    {
                        insideDomain = true;

                        previousCell[i] = currentCell[i];
                        currentCell[i] = neighbourCellNo;

                        break;
                    }
                    if(current_depth + 1 < searchDepth)
                        // add the neighbour cell to the queue
                        cell_queue.push({neighbourCellNo, current_depth+1});
                }
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
                    isParticleDeposited[i] = 1;
                    continue;
                }
                else
                {
                    // Check for Boundary of the neighbours for particles.
                    cout << "Error" << endl;
                    isParticleDeposited[i] = 1;
                    m_ErrorParticlesCount++;
                    isErrorParticle[i] = 1;
                    continue;
                }
            }

            // Its a Boundary Cell and the cornerID is 21
            // Mark it as escaped. 
            // Do not perform the intersection logic for a particle , which is at verge of escaping
            if (isBoundaryCell && cornerID == 21 ) 
            {
                isParticleDeposited[i] = 1;
                m_EscapedParticlesCount++;
                isEscapedParticle[i] = 1;
                continue;
            }

             // Check if the particle is from a corner shared by two bdids 2 and 0
            if (isBoundaryCell && cornerID == 20)
            {
                #pragma omp critical
                jointID  = m_jointidOfBoundCells[cellNo];  
            }

            // If a paricle escaped from a cell, in which there is only inlet surface then its an error particle
            // This will only happen if there is a backflow
            if(cornerID == 0)
            {
                isParticleDeposited[i] = 1;
                m_ErrorParticlesCount++;
                isErrorParticle[i] = 1;
                continue;
            }
            else // Corner ID is 1 or 2
            {   
                #pragma omp critical
                jointID = m_jointidOfBoundCells[cellNo];

            }
            
            TJoint *Joint = cell->GetJoint(jointID);
            double x1, x2, x3, y1, y2, y3, z1, z2, z3;

            // Get the coordinates of the joint
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
                isParticleDeposited[i] = 1;

                

                position_X[i] = temp2.x;
                position_Y[i] = temp2.y;
                position_Z[i] = temp2.z;

                // if bdId/corner ID is 1, then mark the particle as escaped
                if(cornerID == 1)
                {
                    isParticleDeposited[i] = 1;
                    m_EscapedParticlesCount++;
                    isEscapedParticle[i] = 1;   // Mark the particle as escaped
                    continue;
                }
            }
            else
            {
                // if bdId/corner ID is 1, then mark the particle as escaped
                if(cornerID == 1)
                {
                    isParticleDeposited[i] = 1;
                    m_EscapedParticlesCount++;
                    isEscapedParticle[i] = 1;   // Mark the particle as escaped
                    continue;
                }

                
                // Mark the Particle as deposited
                isParticleDeposited[i] = 1;
                // Ghost particle, Make the vertex as deposition
                position_X[i] = x1;
                position_Y[i] = y1;
                position_Z[i] = z1;
                m_ghostParticlesCount++;

                
            }
        }
        // check if the particle is in any border cell
    }

    // finih the timer
    auto finish = std::chrono::high_resolution_clock::now();

    // calculate the time taken from start to finish
    std::chrono::duration<double> elapsed = finish - start;

    // print the time taken
    std::cout << "[INFORMATION] Time taken for Interpolation: " << elapsed.count() << " s\n";

    double x_difference_position =  diff_dot_product(position_X,position_X_old);
    double y_difference_position =  diff_dot_product(position_Y,position_Y_old);
    double z_difference_position =  diff_dot_product(position_Z,position_Z_old);
    cout << " Difference in position X : " << x_difference_position<< "\n";
    cout << " Difference in position Y : " << y_difference_position<< "\n";
    cout << " Difference in position Z : " << z_difference_position<< "\n";

    int depositedCount = std::count(isParticleDeposited.begin(), isParticleDeposited.end(), 1);
    m_ErrorParticlesCount = std::count(isErrorParticle.begin(), isErrorParticle.end(), 1);
    m_ghostParticlesCount = std::count(isGhostParticle.begin(), isGhostParticle.end(), 1);
    m_StagnantParticlesCount = std::count(isStagnantParticle.begin(), isStagnantParticle.end(), 1);

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
#endif
