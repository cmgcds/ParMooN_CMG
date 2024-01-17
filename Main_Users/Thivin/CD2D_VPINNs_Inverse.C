// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <MainUtilities.h>

#include <stdlib.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <sstream>
#include <string>

#ifdef _TEST                               // Do not change ..!!!!
#include "../Examples/CD_2D/SineLaplace.h" // smooth sol in unitsquares
#else
// #include "../Examples/CD_2D/SineLaplace.h"
// #include "../Main_Users/Thivin/fastVPINNs/gear.h"
#include "../Main_Users/Thivin/fastVPINNs/singularly_perturbled.h"
#endif

// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/CD_2D/Hemker1996.h" // circle in a channel

// #include "../Examples/CD_2D/TwoInteriorLayers.h" // smooth sol in unitsquares
// #include "../Examples/CD_2D/furnace.h"
// #include "../Main_Users/Sashi/TCD_2D/Hemker.h"

// =======================================================================
// main program
// =======================================================================
int main(int argc, char *argv[])
{
    //  declaration of database, you need this in every program
    int i, ORDER, N_Cells, N_DOF, img = 0, Disctype;

    double *sol, *rhs, t1, t2, errors[4];

    char *VtkBaseName;

    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();
    TCollection *coll;
    TDomain *Domain;
    TFESpace2D *Scalar_FeSpace, *fesp[1];
    TFEFunction2D *Scalar_FeFunction;
    TSystemCD2D *SystemMatrix;
    TOutput2D *Output;
    TAuxParam2D *aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

    std::ostringstream os;
    os << " ";

    // set variables' value in TDatabase using argv[1] (*.dat file)
    Domain = new TDomain(argv[1]);

    // set PROBLEM_TYPE to CD if not yet set
    if (TDatabase::ParamDB->PROBLEM_TYPE == 0)
        TDatabase::ParamDB->PROBLEM_TYPE = 1;
    // open OUTFILE, this is where all output is written to (addionally to console)
    OpenFiles();

    // write all Parameters to the OUTFILE (not to console) for later reference
    Database->WriteParamDB(argv[0]);
    ExampleFile();

    /* include the mesh from a mesh generator, for a standard mesh use the
     * build-in function. The GEOFILE describes the boundary of the domain. */
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
        OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
        OutPut("GMSH used for meshing !!!" << endl);

    }                                            // gmsh mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 2) // triangle mesh
    {
        OutPut("Triangle.h used for meshing !!!" << endl);
        TriaReMeshGen(Domain);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }

#if defined(__HEMKER__) || defined(__BEAM__)
    TriaReMeshGen(Domain);
    TDatabase::ParamDB->UNIFORM_STEPS = 0;
#endif

    // refine grid up to the coarsest level

    for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    // write grid into an Postscript file
    if (TDatabase::ParamDB->WRITE_PS)
        Domain->PS("Domain.ps", It_Finest, 0);

    // create output directory, if not already existing
    if (TDatabase::ParamDB->WRITE_VTK)
        mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

    //=========================================================================
    // construct all finite element spaces
    //=========================================================================
    cout << "__THIVIN__ ANSATZ_ORDER : " << TDatabase::ParamDB->ANSATZ_ORDER << endl;
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

    // a collection is basically only an array of cells, which is needed to create
    // a finite element space
    coll = Domain->GetCollection(It_Finest, 0);
    // print out some information about the mesh
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells : " << N_Cells << endl);

    // create fespace for scalar equation
    Scalar_FeSpace = new TFESpace2D(coll, (char *)"name", (char *)"description", BoundCondition, ORDER, NULL);
    // print out some information on the finite element space
    N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    OutPut("dof all      : " << setw(10) << N_DOF << endl);

    //======================================================================
    // construct all finite element functions
    //======================================================================
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    // set solution and right hand side vectors to zero
    memset(sol, 0, N_DOF * SizeOfDouble);
    memset(rhs, 0, N_DOF * SizeOfDouble);

    double* eps = new double[N_DOF];
    // create a finite element function
    Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char *)"C", (char *)"C", sol, N_DOF);

    TFEFunction2D* eps_feFunction = new TFEFunction2D(Scalar_FeSpace, (char *)"eps", (char *)"eps", eps, N_DOF);


    // from fespace get all the DOF Coordinates
    double* x_cord = new double[N_DOF];
    double* y_cord = new double[N_DOF];
    Scalar_FeSpace->GetDOFPosition(x_cord, y_cord);

    // Loop over all DOF and obtain the corresponding eps values from the bilinearcoeffs
    for (int i = 0; i < N_DOF; i++)
    {
        double x = x_cord[i];
        double y = y_cord[i];
        
        int n_param = 1;
        // create a double pointer for param, where the first dimension is n_param and inner dimension is 4
        double** param = new double*[n_param];
        for (int j = 0; j < n_param; j++)
        {
            param[j] = new double[6];
        }

        BilinearCoeffs(n_param, &x, &y, NULL, param);

        eps[i] = param[0][0];
    }

    //======================================================================
    // SystemMatrix construction and solution
    //======================================================================
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    Disctype = TDatabase::ParamDB->DISCTYPE;
    SystemMatrix = new TSystemCD2D(Scalar_FeSpace, Disctype, DIRECT);

    // initilize the system matrix with the functions defined in the example
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);

    // assemble the system matrix with given aux, sol and rhs
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass with NULL
    SystemMatrix->Assemble(NULL, sol, rhs);

    // Solve the system
    t1 = GetTime();
    SystemMatrix->Solve(sol, rhs);
    t2 = GetTime();
    OutPut("time for solving: " << t2 - t1 << endl);

    // Generate a uniformly spaced 2D Grid with mesh size as 0.01 from 0 t 1

    


    // // std::string vpinnfilename = "vpinn_test_points_circle.csv";
    // std::string vpinnfilename = "test_points_gear.csv";

    // // Open this csv and read the 2d points seperated by a space to a vector of doubles
    // std::ifstream vpinnfile(vpinnfilename);

    // // if file is not found, throw an error
    // if (!vpinnfile)
    // {
    //     std::cout << "Error opening file " << vpinnfilename << std::endl;
    //     exit(1);
    // }

    // std::vector<double> vpinnpoints_x;
    // std::vector<double> vpinnpoints_y;
    // std::string line;

    // while (std::getline(vpinnfile, line))
    // {
    //     std::istringstream iss(line);
    //     double x, y;
    //     if (!(iss >> x >> y))
    //     {
    //         break;
    //     }
    //     vpinnpoints_x.push_back(x);
    //     vpinnpoints_y.push_back(y);
    // }




    // std::string vpinnoutputfilename = "fem_output_gear_forward_sin.csv";
    // std::ofstream vpinnoutputfile(vpinnoutputfilename);

    // // Loop over all the points and find the corresponding solution and eps values
    // for (int i = 0; i < vpinnpoints_x.size(); i++)
    // {
    //     if(i % 10000 == 0)
    //     {
    //         std::cout << "Processing point " << i << std::endl;
    //     }
    //     double x = vpinnpoints_x[i];
    //     double y = vpinnpoints_y[i];

    //     int n_param = 1;
    //     // create a double pointer for param, where the first dimension is n_param and inner dimension is 4
    //     double** param = new double*[n_param];
    //     for (int j = 0; j < n_param; j++)
    //     {
    //         param[j] = new double[6];
    //     }

    //     BilinearCoeffs(n_param, &x, &y, NULL, param);

    //     double eps = param[0][0];

    //     double values[4];
    //     Scalar_FeFunction->FindGradient(x, y, values);

    //     // Write the output to the file
    //     vpinnoutputfile << x << "," << y << "," << values[0] << "," << values[1] << "," << values[2] << "," << eps << endl;
        
    // }

    // std::cout << "Output written to " << vpinnoutputfilename << std::endl;

    // create a file to store the vpinn output
    // Assume the limits and n_test_points are given
    double x_start = 0., x_end = 1.0;
    double y_start = 0., y_end = 1.0;

    int n_test = (int)(TDatabase::ParamDB->P15);

    int n_test_points_x = n_test;
    int n_test_points_y = n_test;

    // Compute the step sizes
    double x_step = (x_end - x_start) / (n_test_points_x - 1);
    double y_step = (y_end - y_start) / (n_test_points_y - 1);

    // Create the x and y vectors
    std::vector<double> x(n_test_points_x), y(n_test_points_y);
    for (int i = 0; i < n_test_points_x; ++i) {
        x[i] = x_start + i * x_step;
    }
    for (int i = 0; i < n_test_points_y; ++i) {
        y[i] = y_start + i * y_step;
    }

    // Create the meshgrid and stack the points
    std::vector<double> test_points_x(n_test_points_x * n_test_points_y);
    std::vector<double> test_points_y(n_test_points_x * n_test_points_y);
    for (int i = 0; i < n_test_points_x; ++i) {
        for (int j = 0; j < n_test_points_y; ++j) {
            test_points_x[i * n_test_points_y + j] = x[i];
            test_points_y[i * n_test_points_y + j] = y[j];
        }
    }

    std::string filename_prefix ;

    if (TDatabase::ParamDB->DISCTYPE == 1)
    {
        filename_prefix = "fem_output_eps_";
    }
    else
    {
        filename_prefix = "supg_fem_output_eps_";
    }


    double eps_val = TDatabase::ParamDB->PE_NR;
    std::ostringstream streamObj;
    streamObj << eps_val;
    std::string strEps = streamObj.str();

    std::string vpinn_output_file_name = filename_prefix  +  strEps + ".csv";

    
    std::ofstream vpinnoutputfile(vpinn_output_file_name);

    // Loop over all the points and find the corresponding solution and eps values
    for (int i = 0 ; i < test_points_x.size(); i++)
    {
        if(i % 10000 == 0)
        {
            std::cout << "Processing point " << i << std::endl;
        }
        double x = test_points_x[i];
        double y = test_points_y[i];

        int n_param = 1;
        // create a double pointer for param, where the first dimension is n_param and inner dimension is 4
        double** param = new double*[n_param];
        for (int j = 0; j < n_param; j++)
        {
            param[j] = new double[6];
        }

        BilinearCoeffs(n_param, &x, &y, NULL, param);

        double eps = param[0][0];

        double values[4];
        Scalar_FeFunction->FindGradient(x, y, values);
        vpinnoutputfile << x << "," << y << "," << values[0] << "," << values[1] << "," << values[2] << "," << eps << endl;

        
    }

    std::cout << "Output written to " << vpinn_output_file_name << std::endl;

    //======================================================================
    // produce outout
    //======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    Output = new TOutput2D(1, 1, 0, 0, Domain);
    Output->AddFEFunction(Scalar_FeFunction);
    Output->AddFEFunction(eps_feFunction);

    //     Scalar_FeFunction->Interpolate(Exact);
    if (TDatabase::ParamDB->WRITE_VTK)
    {
        os.seekp(std::ios::beg);
        os << "VTK/" << VtkBaseName << "_"<< strEps << "." << img << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        // if (img < 10)
        //     os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
        // else if (img < 100)
        //     os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
        // else if (img < 1000)
        //     os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
        // else if (img < 10000)
        //     os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
        // else
        //     os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;
    }

    //======================================================================
    // measure errors to known solution
    // If an exact solution is not known, it is usually set to be zero, so that
    // in such a case here only integrals of the solution are computed.
    //======================================================================
    if (TDatabase::ParamDB->MEASURE_ERRORS)
    {
        fesp[0] = Scalar_FeSpace;
        aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

        Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                                     BilinearCoeffs, aux, 1, fesp, errors);

        delete aux;

        OutPut("L2: " << errors[0] << endl);
        OutPut("H1-semi: " << errors[1] << endl);
        OutPut("SD: " << errors[2] << endl);

    } // if(TDatabase::ParamDB->MEASURE_ERRORS)

#ifdef _TEST
    cout << "TESTVALUE_CD2D " << errors[0] << endl;
#endif

    delete[] sol;
    delete[] rhs;
    delete coll;
    CloseFiles();
    return 0;
} // end main
