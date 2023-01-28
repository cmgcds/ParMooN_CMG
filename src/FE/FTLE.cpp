#include <AllClasses.h>
#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <DeformMesh2D.h>
#include <LinAlg.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <FEVectFunct2D.h>
#include <Constants.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <NodalFunctional2D.h>
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

#include<omp.h>

// PRogress Bar

// COnstructor
FTLE::FTLE(TFESpace2D *ftleFespace, TFEVectFunct2D *FEVectFunction, int searchDepth, TFEVectFunct2D *uCombined, TFEVectFunct2D *vCombined)
{

    fespace = ftleFespace;
    VelocityVectFunction = FEVectFunction;
    N_Particles = ftleFespace->GetN_DegreesOfFreedom();

    TCollection *coll = ftleFespace->GetCollection();
    N_cells = fespace->GetN_Cells();
    int N_Points;
    double *xi, *eta;
    double X[MaxN_BaseFunctions2D], Y[MaxN_BaseFunctions2D];
    RefTrans2D RefTrans, *RefTransArray;

    u_Vectfunc_Combined = uCombined;
    v_Vectfunc_Combined = vCombined;

    // Construct the arrays for position and the cellId's
    x_pos = new double[N_Particles]();
    y_pos = new double[N_Particles]();
    cellIds = new int[N_Particles]();
    x_pos_initial = new double[N_Particles]();
    y_pos_initial = new double[N_Particles]();
    isParticleOnBoundary = new int[N_Particles]();

    for (int i = 0; i < N_Particles; i++)
        isParticleOnBoundary[i] = 0;

    Cell2DOFMapping.resize(N_cells, std::vector<int>());
    DOF2CellMapping.resize(N_Particles, std::vector<int>());
    DOF_Neibhours.resize(N_cells, std::vector<int>());
    DOF_Neibhours_depth.resize(N_cells, std::vector<int>());

    // INitialise particles inside the domain to be true.
    particleInsideDomain.resize(N_Particles, true);

    // resize FTLE
    FTLEValues.resize(N_Particles, 0.0);

    neibhourLeft.resize(N_Particles);
    neibhourRight.resize(N_Particles);
    neibhourTop.resize(N_Particles);
    neibhourBottom.resize(N_Particles);

    // Updates the position of All particles in the domain.
    fespace->GetDOFPosition(x_pos, y_pos, cellIds);
    fespace->GetDOFPosition(x_pos_initial, y_pos_initial, cellIds);

    // Get the cell Information
    int *GlobalNumbers = fespace->GetGlobalNumbers();
    int *BeginIndex = fespace->GetBeginIndex();
    const int *TmpLen;
    const int *TmpEV;
    int MaxLen;
    BoundCond Bdcond;
    TBoundEdge *Bdface; // Pointer to Boundary Face in 3D Cell
    TBoundComp *BoundComp;
    TVertex *currentVertex;

    for (int cellId = 0; cellId < N_cells; cellId++)
    {
        TBaseCell *currentCell = coll->GetCell(cellId);
        int *GlobalDOF = GlobalNumbers + BeginIndex[cellId];
        FE2D elementId = fespace->GetFE2D(cellId, currentCell);
        TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
        TFEDesc2D *fedesc = element->GetFEDesc2D();
        // currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        int N_Joints = currentCell->GetN_Joints();

        for (int jointId = 0; jointId < N_Joints; jointId++)
        {
            TJoint *Joint = currentCell->GetJoint(jointId);
            if (Joint->GetType() == BoundaryEdge)
            {
                int *JointDOF = fedesc->GetJointDOF(jointId);
                for (int vert = 0; vert < fedesc->GetN_JointDOF(); vert++)
                {
                    // int local_vertex_no     =   TmpFV[jointId*MaxLen + vert];
                    int glob_vertex_no = GlobalDOF[JointDOF[vert]];
                    isParticleOnBoundary[glob_vertex_no] = 1;
                    // Update the position of the DOF based on the cell structure.
                }
            }
        }
    }

    // for ( int i=0 ; i < N_Particles; i++)
    // {
    //     cout << " DOF : "<<setw(3) << i <<setw(3) << " x : " <<setw(8)<< x_pos_initial[i] << " y: "<<setw(8)<< y_pos_initial[i] << " Cell : " <<setw(5)<< cellIds[i] << " Boundary : " <<setw(8)<< isParticleOnBoundary[i]<<endl;
    // }

    int *GlobalDOFArray = fespace->GetGlobalNumbers();
    int *BeginIndexArray = fespace->GetBeginIndex();

    // Get Reftrans Array
    RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

    cout << " N_cells : " << N_cells << endl;
    cout << " N_Particles : " << N_Particles << endl;

    // Create a vector of vector with size N_Particles;

    for (int cellno = 0; cellno < N_cells; cellno++)
    {

        int *globalDOFCell = fespace->GetGlobalDOF(cellno);
        // Get the cell , Reference transformation and its basis functions to get the local DOF>
        TBaseCell *cell = coll->GetCell(cellno);
        FE2D FEid = fespace->GetFE2D(cellno, cell);
        RefTrans = RefTransArray[FEid];
        BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEid);

        TNodalFunctional2D *nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(FEid);
        nf->GetPointsForAll(N_Points, xi, eta);

        // cout << " N_Points : "<< N_Points <<endl;

        // printf(" hello");
        for (int j = 0; j < N_Points; j++)
        {
            int dof = globalDOFCell[j];
            Cell2DOFMapping[cellno].push_back(dof);
            // printf("%d ", j);
            DOF2CellMapping[dof].push_back(cellno);
            // printf("%d ", j);
        }
    }

    // Create DOF Neibhous
    for (int cellNo = 0; cellNo < N_cells; cellNo++)
    {
        std::set<int> ss;
        std::vector<int> v;

        for (int j = 0; j < Cell2DOFMapping[cellNo].size(); j++)
        {
            int dof = Cell2DOFMapping[cellNo][j];
            // cout << " outerLoop : DOF no : "<< dof <<endl;
            for (int cellofDof = 0; cellofDof < DOF2CellMapping[dof].size(); cellofDof++)
            {
                // cout << DOF2CellMapping[dof][cellofDof] << " ";
                int NeibhcellNo = DOF2CellMapping[dof][cellofDof];
                if (NeibhcellNo != cellNo)
                    ss.insert(DOF2CellMapping[dof][cellofDof]);
            }
            // cout<<endl;
        }
        v.assign(ss.begin(), ss.end());
        DOF_Neibhours[cellNo] = v;
        ss.clear();
    }

    // cout << "MaxValue : " << *std::max_element(cellIds,cellIds+N_Particles)     <<endl;

    for (int cellId = 0; cellId < N_cells; cellId++)
    {
        std::vector<int> immediateNeibhCells = DOF_Neibhours[cellId];
        int start = 0;
        int prevSize = immediateNeibhCells.size();
        // cout << "------- Cell : " << cellId << endl;
        for (int depth = 0; depth < searchDepth; depth++)
        {

            for (int j = start; j < prevSize; j++)
            {
                // Get neibhours of these cells
                int neibhCell = immediateNeibhCells[j];

                for (int k = 0; k < DOF_Neibhours[neibhCell].size(); k++)
                {
                    int cell = DOF_Neibhours[neibhCell][k];
                    if (cell != cellId && (std::find(immediateNeibhCells.begin(), immediateNeibhCells.end(), cell) == immediateNeibhCells.end()))
                        immediateNeibhCells.push_back(cell);
                }

                // for(int i=0 ; i < immediateNeibhCells.size(); i++)
                // cout << " " << immediateNeibhCells[i] ;
                // cout <<endl;
            }
            start = prevSize;
            prevSize = immediateNeibhCells.size();
        }

        DOF_Neibhours_depth[cellId] = immediateNeibhCells;
    }

    // for( int cellNo = 0; cellNo < N_cells; cellNo++ )
    // {
    //     cout << "Cell mo : "<< cellNo << " --> SIZE : " << DOF_Neibhours_depth[cellNo].size() << " \n";
    //     for( int dof = 0 ; dof < DOF_Neibhours_depth[cellNo].size() ; dof++ )
    //     {
    //         cout << "  " << DOF_Neibhours_depth[cellNo][dof];
    //     }
    //     cout<<endl;
    // }

    // exit(0);

    cout << FGRN(" [Sucess] Cell Neibhours based on DOF - Successfully Generated ") << endl;

    // Calculate the neibhours // for(int particleId = 0 ; particleId < 3; particleId++)
    for (int particleId = 0; particleId < N_Particles; particleId++)
    {
        int currCell = cellIds[particleId];
        // Get all neibhours of the DOF ( )
        std::vector<int> NeibDOFs;

        // cout << " CurrCell : " << currCell <<endl;

        for (int k = 0; k < Cell2DOFMapping[currCell].size(); k++)
        {
            int dof = Cell2DOFMapping[currCell][k];
            // Particle is not same as the current particle
            if (!(abs(x_pos[particleId] - x_pos[dof]) < 1e-9 && (abs(y_pos[particleId] - y_pos[dof])) < 1e-9))
            {
                NeibDOFs.push_back(dof);
            }
        }

        for (int j = 0; j < DOF_Neibhours[currCell].size(); j++)
        {
            int cell = DOF_Neibhours[currCell][j]; // Gets cell no of neibhouring cells
            for (int k = 0; k < Cell2DOFMapping[cell].size(); k++)
            {
                int dof = Cell2DOFMapping[cell][k];
                if (!(abs(x_pos[particleId] - x_pos[dof]) < 1e-9 && (abs(y_pos[particleId] - y_pos[dof])) < 1e-9))
                    NeibDOFs.push_back(dof);
            }
        }

        // cout << " Collection of All nebhours " <<endl;
        // std::for_each(NeibDOFs.begin(),NeibDOFs.end(),[&](auto i){cout << "DOF : " << i << " xpos: " << x_pos[i] << " , " << y_pos[i] <<endl;});
        // cout<<endl;

        std::unordered_set<int> s1(NeibDOFs.begin(), NeibDOFs.end());
        NeibDOFs.assign(s1.begin(), s1.end());

        double minLeft = 10000;
        double minBott = 10000;
        const double par_x = x_pos[particleId];
        const double par_y = y_pos[particleId];
        double samelineTolerance = 0.001;
        int leftDOF;
        bool leftFound = false;
        // cout << " PARTIcle : " << par_x << " , " << par_y <<endl;
        // Check for left neibhours
        for (int i = 0; i < NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x = x_pos[NeibDOFs[i]];
            double y = y_pos[NeibDOFs[i]];
            double left = par_x - x;
            double bott = abs(par_y - y);

            // cout << " ( " <<x <<","<<y<<") "<< " -- left : "<< left << " bott : "<<bott << " minBottom : " << bott  << "  min left : " << minLeft<<endl;

            if ((left > 0) && (left <= minLeft) && (abs(bott) <= minBott))
            {
                //  cout << " &&&&&&&( " <<x <<","<<y<<") "<<endl;
                leftDOF = dof;
                minLeft = left;
                minBott = bott;
                leftFound = true;
            }

            if (leftFound)
                neibhourLeft[particleId] = leftDOF;
            else
                neibhourLeft[particleId] = -9999; // LEFT CORNER DOF
        }
        // cout << " Neibhour Left : " << neibhourLeft[particleId]  << " x: " << x_pos[neibhourLeft[particleId]] << " y : " << y_pos[neibhourLeft[particleId]]<<endl;

        // FIND RIGHT DOF
        double minRight = 100000;
        minBott = 10000;

        samelineTolerance = 0.001;
        int rightDOF;
        bool rightFound = false;
        // Check for left neibhours
        for (int i = 0; i < NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x = x_pos[NeibDOFs[i]];
            double y = y_pos[NeibDOFs[i]];
            double right = x - par_x;
            double bott = abs(par_y - y);

            if ((right > 0) && (right <= minRight) && (abs(bott) <= minBott))
            {
                rightDOF = dof;
                minRight = right;
                minBott = bott;
                rightFound = true;
            }

            if (rightFound)
                neibhourRight[particleId] = rightDOF;
            else
                neibhourRight[particleId] = -9999; // LEFT CORNER DOF
        }
        // cout << " Neibhour Right : " << neibhourRight[particleId]  << " x: " << x_pos[neibhourRight[particleId]] << " y : " << y_pos[neibhourRight[particleId]]<<endl;

        // FIND TOP DOF
        double minTop = 100000;
        double minSide = 10000;
        samelineTolerance = 0.001;
        int topDOF;
        bool topFound = false;
        // Check for left neibhours
        for (int i = 0; i < NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x = x_pos[NeibDOFs[i]];
            double y = y_pos[NeibDOFs[i]];
            double top = y - par_y;
            double side = abs(par_x - x);

            if ((top > 0) && (top <= minTop) && (abs(side) <= minSide))
            {
                topDOF = dof;
                minTop = top;
                minSide = side;
                topFound = true;
            }

            if (topFound)
                neibhourTop[particleId] = topDOF;
            else
                neibhourTop[particleId] = -9999; // LEFT CORNER DOF
        }
        // cout << " Neibhour Top : " << neibhourTop[particleId]  << " x: " << x_pos[neibhourTop[particleId]] << " y : " << y_pos[neibhourTop[particleId]]<<endl;

        // FIND BOTTOM DOF
        double minBottom = 100000;
        minSide = 10000;

        samelineTolerance = 0.001;
        int bottomDOF;
        bool bottomFound = false;
        // Check for left neibhours
        for (int i = 0; i < NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x = x_pos[NeibDOFs[i]];
            double y = y_pos[NeibDOFs[i]];
            double bottom = par_y - y;
            double side = abs(par_x - x);

            if ((bottom > 0) && (bottom <= minBottom) && (abs(side) < minSide))
            {
                bottomDOF = dof;
                minBottom = bottom;
                minSide = side;
                bottomFound = true;
            }

            if (bottomFound)
                neibhourBottom[particleId] = bottomDOF;
            else
                neibhourBottom[particleId] = -9999; // LEFT CORNER DOF
        }

        // cout << " Neibhour Top : " << neibhourBottom[particleId]  << " x: " << x_pos[neibhourBottom[particleId]] << " y : " << y_pos[neibhourBottom[particleId]]<<endl
    }

    cout << FGRN(" [Sucess] Directional Neibhours Obtained Successfully ") << endl;
}

double *FTLE::computeFTLE(double timeStep, int T, int startT)
{

    double *xpos_final = new double[N_Particles]();
    double *ypos_final = new double[N_Particles]();
    int *cellTracker = new int[N_Particles]();

    // Copy the intial position to the given particles.
    for (int i = 0; i < N_Particles; i++)
    {
        xpos_final[i] = x_pos_initial[i];
        ypos_final[i] = y_pos_initial[i];
        cellTracker[i] = cellIds[i];
    }

    for (int i = 0; i < N_Particles; i++)
    {
        particleInsideDomain[i] = true;
    }

    //   for ( int i=0 ; i < N_Particles; i++)
    //     cout << " DOF : "<< i << " x : " << xpos_final[i] << " y: "<< ypos_final[i] << " Cell : " << cellTracker[i]<<endl;


    //  Set all FTLE values to be zero
    std::fill(FTLEValues.begin(), FTLEValues.end(), 0.0);

    // std::string filename = "position_" + std::to_string(0) + ".csv" ;

    // std::ofstream myfile(filename);
    //     for (int i = 0 ; i < FTLEValues.size();i++)
    // {
    //     myfile << xpos_final[i] <<","<<ypos_final[i] <<"," << FTLEValues[i] <<endl;
    // }

    // myfile.close();

    cout << " Start TIme : " << startT << " End : " << startT + T - 1 << endl;
    double t1 = omp_get_wtime();
    for (int time = startT; time < startT + T - 1; time++)
    {
        int outside = 0;
        int otherCountTemp = 0;
        
        TFEFunction2D* uVelocity = u_Vectfunc_Combined->GetComponent(time);
        TFEFunction2D* vVelocity = v_Vectfunc_Combined->GetComponent(time);
        TFEFunction2D* uVelocity_new = u_Vectfunc_Combined->GetComponent(time+1);
        TFEFunction2D* vVelocity_new = v_Vectfunc_Combined->GetComponent(time+1);

        

        
        omp_set_num_threads(24);
        #pragma omp parallel for schedule(static,8) shared(particleInsideDomain, fespace, u_Vectfunc_Combined, v_Vectfunc_Combined, \
                                        xpos_final, ypos_final, cellTracker, N_cells, timeStep,DOF_Neibhours_depth,uVelocity,vVelocity,uVelocity_new,vVelocity_new)
        for (int i = 0; i < N_Particles; i++)
        {
            if (!particleInsideDomain[i])
            {
                continue; // If particle is not inside the domain , then proceed to next particle
                otherCountTemp++;
            }

            int currentCellNo = cellTracker[i];

            std::vector<int> neibhours ;
            // neibhours.push_back(currentCellNo);
            for (int k = 0; k < N_cells; k++)
                neibhours.push_back(k);

            // neibhours.push_back(currentCellNo);
            bool insideDomain = false;
            TBaseCell *cell;
            for (int cellId = 0; cellId < neibhours.size(); cellId++)
            {
                cell = fespace->GetCollection()->GetCell(neibhours[cellId]);

                if (cell->PointInCell(xpos_final[i], ypos_final[i]))
                {
                    insideDomain = true;
                    cellTracker[i] = cellId;
                    break;
                }
            }

            if (!insideDomain)
            {
                particleInsideDomain[i] = false;
                outside++;
                continue;
            }

            currentCellNo = cellTracker[i];

            double uVal[3];
            double vVal[3];

            double uValNew[3];
            double vValNew[3];

            double uValCurrent, vValCurrent;
            double uValNext = 0, vValNext = 0;

            // u_Vectfunc_Combined->GetComponent(time)->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], uVal);
            // v_Vectfunc_Combined->GetComponent(time)->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], vVal);
            uVelocity->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], uVal);
            vVelocity->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], vVal);
            // uVal[0] = 0;
            // vVal[0] = 0;
            uValCurrent = uVal[0];
            vValCurrent = vVal[0];
            
            // u_Vectfunc_Combined->GetComponent(time+1)->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], uValNew);
            // v_Vectfunc_Combined->GetComponent(time+1)->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], vValNew);
            uVelocity_new->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], uValNew);
            vVelocity_new->FindValueLocal_Parallel(cell, currentCellNo, xpos_final[i], ypos_final[i], vValNew);
            // uValNew[0] = 0;
            // vValNew[0] = 0;
            uValNext = uValNew[0];
            vValNext = vValNew[0];


            // cout << 0.5 * (uValCurrent + uValNext) << ", " << uValCurrent << ", "<< vValCurrent << " , " << uValNext << ", " << vValNext <<endl;
            xpos_final[i] += timeStep * 0.5 * (uValCurrent + uValNext);
            ypos_final[i] += timeStep * 0.5 * (vValCurrent + vValNext);

            neibhours.clear();
 
        }

    }

    cout << " time elapsed ; " << omp_get_wtime() - t1 <<endl;

    t1 = omp_get_wtime();

    int counttt = 0;
    int outside = 0;

    // #pragma omp parallel for schedule(static,8) shared(particleInsideDomain, neibhourLeft, neibhourRight, neibhourTop, neibhourBottom, \
    //                             xpos_final, ypos_final, FTLEValues,outside)
    // After computing Displacement, Compute the FTLE Coefficients
    for (int i = 0; i < N_Particles; i++)
    {
        if (!particleInsideDomain[i] )
        {
            // #pragma omp critical
            outside++;
            continue; // If particle is not inside the domain , then proceed to next particle
        }
        int n_left = neibhourLeft[i];
        int n_right = neibhourRight[i];
        int n_top = neibhourTop[i];
        int n_bottom = neibhourBottom[i];

        double a11 = 0, a12 = 0, a21 = 0, a22 = 0;
        // Do not compute the FTLE for Border nodes.
        if (n_left != -9999 && n_right != -9999)
        // if (!isParticleOnBoundary[i] )
        {
            a11 = (xpos_final[n_right] - xpos_final[n_left]) / (x_pos_initial[n_right] - x_pos_initial[n_left]);
            a21 = (ypos_final[n_right] - ypos_final[n_left]) / (x_pos_initial[n_right] - x_pos_initial[n_left]);

            counttt++;
            
        }

        if (n_top != -9999 && n_bottom != -9999)
        // if (!isParticleOnBoundary[i])
        {
            a12 = (xpos_final[n_top] - xpos_final[n_bottom]) / (y_pos_initial[n_top] - y_pos_initial[n_bottom]);
            a22 = (ypos_final[n_top] - ypos_final[n_bottom]) / (y_pos_initial[n_top] - y_pos_initial[n_bottom]);
        }

        double *C = new double[2 * 2]();
        // Matrix Matrix Multiplication
        // C = AT*A

        C[0] = a11 * a11 + a21 * a21;
        C[1] = a11 * a12 + a21 * a22;
        C[2] = a11 * a12 + a21 * a22;
        C[3] = a21 * a21 + a22 * a22;

        // Compute SVD
        int m = 2, n = 2;
        int lda = m;
        int ldu = m;
        int ldvt = m;
        int info;
        int lwork = -1;
        double wkopt;
        double *work;

        double a[] = {C[0], C[2], C[1], C[3]}; // Deformation Tensor value in column major
        double s[n], u[ldu * m], vt[ldvt * n];

        dgesvd("All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);
        lwork = (int)wkopt;
        work = (double *)malloc(lwork * sizeof(double));
        /* Compute SVD */
        dgesvd("All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,
               &info);

        if (info > 0)
        {
            printf("The algorithm computing SVD failed to converge.\n");
            exit(1);
        }

        delete[] C;
        // printf("Pos vALUES : %f, %f , %f, %f\n",x_pos_initial[i],y_pos_initial[i],xpos_final[i],ypos_final[i]);
        // printf("mAT vALUES : %f, %f , %f, %f\n",a11,a12,a21,a22);
        // printf("mAT vALUES : %f, %f , %f, %f\n",C[0],C[1],C[2],C[3]);
        // printf(" Sign Values : %f, %f\n",s[0],s[1]);
        double maxSvd = 0.0;
        if (abs(s[0]) > abs(s[1]))
            maxSvd = s[0];
        else
            maxSvd = s[1];

        if (abs(maxSvd - 0.0) > 1e-9)
            // FTLEValues[i] = (1.0 / abs(T * timeStep)) * log(sqrt(maxSvd));
            FTLEValues[i] =  log((maxSvd));
        // FTLEValues[i] = sqrt(maxSvd);
        else
            FTLEValues[i] = 0.0;

        // exit(0);
    }

    cout << " time elapsed ; " << omp_get_wtime() - t1 <<endl;

    std::string filename = "position_" + std::to_string(startT)+ ".csv";

    std::ofstream myfile(filename);
    for (int i = 0; i < FTLEValues.size(); i++)
    {
        myfile << xpos_final[i] << "," << ypos_final[i] << "," << x_pos_initial[i] << "," << y_pos_initial[i] << "," << FTLEValues[i]  << ","
           << xpos_final[neibhourRight[i]] <<"," << ypos_final[neibhourRight[i]] << " , " << xpos_final[neibhourLeft[i]] <<"," << ypos_final[neibhourLeft[i]] << "," 
           << xpos_final[neibhourTop[i]] <<"," << ypos_final[neibhourTop[i]] << "," << xpos_final[neibhourBottom[i]] << "," << ypos_final[neibhourBottom[i]] << ","
            << x_pos_initial[neibhourRight[i]] <<"," << y_pos_initial[neibhourRight[i]] << " , " << x_pos_initial[neibhourLeft[i]] <<"," << y_pos_initial[neibhourLeft[i]] << "," 
           << x_pos_initial[neibhourTop[i]] <<"," << y_pos_initial[neibhourTop[i]] << "," << x_pos_initial[neibhourBottom[i]] << "," << y_pos_initial[neibhourBottom[i]] << endl;
           
    }

    myfile.close();

    cout << " Total : " << N_Particles << " Counted : " << N_Particles - outside << " Outside : " << outside << " Percentage : " << (double)double(counttt) / double(N_Particles) << endl;
    cout << " FTLE Norm : " << sqrt(Ddot(N_Particles, &FTLEValues[0], &FTLEValues[0])) << endl;
    mkdir("FTLE", 0777);
    // std::string name = "FTLE/data_" + std::to_string(int(startT)) + ".csv";
    // write_vtk(name);

    delete [] xpos_final;
    delete [] ypos_final;
    delete [] cellTracker;
}

// Write the FTLE valyes to a VTK, So that they can be visualised using paraview
void FTLE::write_vtk(std::string fileName)
{
    std::ofstream myfile(fileName);
    std::ostringstream os;
    os << " ";

    if (!myfile.is_open())
    {
        cout << FRED(" [Error] File cannot be opened for Writting ") << endl;
        cout << " FileName : FTLE.cpp "
             << " Function : write_vtk() " << endl;
        exit(0);
    }

    for (int i = 0; i < FTLEValues.size(); i++)
    {
        myfile << x_pos[i] << "," << y_pos[i] << "," << FTLEValues[i] << endl;
    }

    myfile.close();
}

// Write the FTLE valyes to a VTK, So that they can be visualised using paraview
// void FTLE::write_tecplotDatFile()
// {

// }