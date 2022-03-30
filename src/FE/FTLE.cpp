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
#include <stdlib.h>
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
#include <set>
#include <unordered_set>
#include <algorithm>
#include<colors.h>


// COnstructor
FTLE::FTLE(TFEVectFunct2D *FEVectFunction, int searchDepth)
{
    fespace = FEVectFunction->GetFESpace2D();
    TCollection *coll = fespace->GetCollection();
    VelocityVectFunction = FEVectFunction;
    int N_DOF = fespace->GetN_DegreesOfFreedom();
    int N_cells = fespace->GetN_Cells();
    int N_Points;
    double *xi, *eta;
    double X[MaxN_BaseFunctions2D], Y[MaxN_BaseFunctions2D];
    RefTrans2D RefTrans, *RefTransArray;

    
    // Construct the arrays for position and the cellId's
    x_pos = new double[N_DOF]();
    y_pos = new double[N_DOF]();
    cellIds = new int[N_DOF]();

    Cell2DOFMapping.resize(N_cells, std::vector<int>());
    DOF2CellMapping.resize(N_DOF, std::vector<int>());
    DOF_Neibhours.resize(N_cells, std::vector<int>());
    DOF_Neibhours_depth.resize(N_cells, std::vector<int>());
    particleInsideDomain.resize(N_DOF,true);
    neibhourLeft.resize(N_DOF);
    neibhourRight.resize(N_DOF);
    neibhourTop.resize(N_DOF);
    neibhourBottom.resize(N_DOF);

    // Updates the position of All particles in the domain.
    fespace->GetDOFPosition(x_pos, y_pos, cellIds);

    int *GlobalDOFArray = fespace->GetGlobalNumbers();
    int *BeginIndexArray = fespace->GetBeginIndex();

    //Get Reftrans Array
    RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

    cout << " N_cells : " << N_cells << endl;
    cout << " N_DOF : " << N_DOF << endl;
    //Create a vector of vector with size N_DOF;

    

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

    // // DOF MAPPING
    // cout << "DOF MAPPING " <<endl;
    // for( int i = 0 ; i < DOF2CellMapping.size() ; i++)
    // {
    //     cout << " DOF : " << i <<" --> ";
    //     for ( int j = 0 ; j < DOF2CellMapping[i].size(); j++)
    //     {
    //         cout <<  ", "<<DOF2CellMapping[i][j] ;
    //     }
    //     cout <<endl;
    // }

    // // CELL MAPPING
    // cout << "CELL MAPPING " <<endl;
    // for( int i = 0 ; i < Cell2DOFMapping.size() ; i++)
    // {
    //     cout << " CELL : " << i <<" --> " ;
    //     for ( int j = 0 ; j < Cell2DOFMapping[i].size(); j++)
    //     {
    //         cout <<  ", "<<Cell2DOFMapping[i][j] ;
    //     }
    //     cout <<endl;
    // }

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

    

    for (int cellId = 0; cellId < N_cells; cellId++)
    {
        std::vector<int> immediateNeibhCells = DOF_Neibhours[cellId];
        int start = 0;
        int prevSize = immediateNeibhCells.size();
        // cout << "------- Cell : " << cellId << endl;
        for (int depth = 0; depth < 2; depth++)
        {

            for (int j = start; j < prevSize; j++)
            {
                // Get neibhours of these cells
                int neibhCell = immediateNeibhCells[j];

                for (int k = 0; k < DOF_Neibhours[neibhCell].size(); k++)
                {
                    int cell = DOF_Neibhours[neibhCell][k];
                    if (cell != cellId && (std::find(immediateNeibhCells.begin(),immediateNeibhCells.end(),cell) == immediateNeibhCells.end() ))
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
    //     cout << "Cell mo : "<< cellNo << " --> ";
    //     for( int dof = 0 ; dof < DOF_Neibhours_depth[cellNo].size() ; dof++ )
    //     {
    //         cout << "  " << DOF_Neibhours_depth[cellNo][dof];
    //     }
    //     cout<<endl;
    // }

    cout<< FGRN(" [Sucess] Cell Neibhours - Successfully Generated ")<<endl;

    // Calculate the neibhours 
    for(int particleId = 0 ; particleId < N_Particles; particleId++)
    {
        int currCell = cellIds[particleId];
        //Get all neibhours of the DOF ( )
        std::vector<int> NeibDOFs;

        for(int k = 0 ; k < Cell2DOFMapping[currCell].size(); k++)     
        {
            int dof = Cell2DOFMapping[currCell][k];
            cout << " Pos : " << x_pos[dof] << ", " <<y_pos[dof]<<endl;
            cout << abs(x_pos[particleId] - x_pos[dof] ) << " , " << (abs(y_pos[particleId] - y_pos[dof]) ) <<endl;
            // Particle is not same as the current particle 
            if(abs(x_pos[particleId] - x_pos[dof] ) > 1e-9 || (abs(y_pos[particleId] - y_pos[dof]) ) > 1e-9)
            {
                NeibDOFs.push_back(dof);
                cout << " PoAdds : " << x_pos[dof] << ", " <<y_pos[dof]<<endl;
                
            }
        }


        for (int j = 0 ; j < DOF_Neibhours[currCell].size(); j++)
        {
            int cell = DOF_Neibhours[currCell][j];  // Gets cell no of neibhouring cells
            for( int k = 0 ; k < Cell2DOFMapping[cell].size(); k++)
            {
                int dof  = Cell2DOFMapping[cell][k];
                if(abs(x_pos[particleId] - x_pos[dof] ) > 1e-9 && (abs(y_pos[particleId] - y_pos[dof] )) > 1e-9)
                    NeibDOFs.push_back(dof);
            }
        }   

        std::unordered_set<int> s1(NeibDOFs.begin(),NeibDOFs.end());
        NeibDOFs.assign(s1.begin(),s1.end());
        
        std::for_each(NeibDOFs.begin(),NeibDOFs.end(),[&](const int a){cout << x_pos[a] <<" " << y_pos[a] <<endl;});
        cout <<endl;

        double minLeft = 100000;
        double minBott  = 10000;
        const double par_x = x_pos[particleId];
        const double par_y = y_pos[particleId];
        double samelineTolerance = 0.001;
        int leftDOF ;
        bool leftFound = false;
        cout << " PARTIcle : " << par_x << " , " << par_y <<endl;
        // Check for left neibhours 
        for( int i = 0 ; i <  NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x  = x_pos[NeibDOFs[i]];
            double y  = x_pos[NeibDOFs[i]];
            double left = par_x - x;
            double bott = abs(par_y - y);
            
            cout << " ( " <<x <<","<<y<<") "<< " -- left : "<< left << " bott : "<<bott <<endl;

            if( abs(x - par_x) > samelineTolerance && (left) && (left < minLeft) && (minBott < bott) )
            {
                leftDOF = dof;
                minLeft = left;
                minBott = bott;
                leftFound = true;
            }

            if(leftFound)
                neibhourLeft[particleId] = leftDOF;
            else
                neibhourLeft[particleId] = -9999;  // LEFT CORNER DOF

        }
        exit(0);

        // FIND RIGHT DOF 
         double minRight = 100000;
         minBott  = 10000;
  
         samelineTolerance = 0.001;
        int rightDOF ;
        bool rightFound = false;
        // Check for left neibhours 
        for( int i = 0 ; i <  NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x  = x_pos[NeibDOFs[i]];
            double y  = x_pos[NeibDOFs[i]];
            double right =  x - par_x;
            double bott = abs(par_y - y);


            if( abs(x - par_x) > samelineTolerance && (right) && (right < minRight) && (minBott < bott) )
            {
                rightDOF = dof;
                minRight = right;
                minBott = bott;
                rightFound = true;
            }

            if(rightFound)
                neibhourRight[particleId] = rightDOF;
            else
                neibhourRight[particleId] = -9999;  // LEFT CORNER DOF
            
        }

        // FIND TOP DOF 
         double minTop = 100000;
        double minSide  = 10000;
        samelineTolerance = 0.001;
        int topDOF ;
        bool topFound = false;
        // Check for left neibhours 
        for( int i = 0 ; i <  NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x  = x_pos[NeibDOFs[i]];
            double y  = x_pos[NeibDOFs[i]];
            double top =  y - par_y;
            double side = abs(par_x - x);


            if( abs(y - par_y) > samelineTolerance && (top) && (top < minTop) && (minSide < side) )
            {
                topDOF = dof;
                minTop = top;
                minSide = side;
                topFound = true;
            }

            if(topFound)
                neibhourRight[particleId] = topDOF;
            else
                neibhourRight[particleId] = -9999;  // LEFT CORNER DOF
            
        }


        // FIND BOTTOM DOF 
         double minBottom = 100000;
        minSide  = 10000;

        samelineTolerance = 0.001;
        int bottomDOF ;
        bool bottomFound = false;
        // Check for left neibhours 
        for( int i = 0 ; i <  NeibDOFs.size(); i++)
        {
            int dof = NeibDOFs[i];
            double x  = x_pos[NeibDOFs[i]];
            double y  = x_pos[NeibDOFs[i]];
            double bottom =   par_y-y;
            double side = abs(par_x - x);


            if( abs(y - par_y) > samelineTolerance && (bottom) && (bottom < minBottom) && (minSide < side) )
            {
                bottomDOF = dof;
                minBottom = bottom;
                minSide = side;
                bottomFound = true;
            }

            if(bottomFound)
                neibhourBottom[particleId] = bottomDOF;
            else
                neibhourBottom[particleId] = -9999;  // LEFT CORNER DOF
            
        }




    }

    for(int i = 0 ; i < N_Particles; i++)
    {
        cout <<  " par : " << i <<" Pos : " << x_pos[i] << " , " << y_pos[i] <<endl;
    }

    // for(int i = 0 ; i < N_Particles; i++)
    // {
    //     cout << "--- Particle " << i << " ------- "<<endl;
    //     cout << " Pos : " << x_pos[i] << " , " << y_pos[i] <<endl;

    //     cout << " Left : " ;
    //     if(neibhourLeft[i])
    //         cout << x_pos[neibhourLeft[i]] << " , " << y_pos[neibhourLeft[i]]  << " DOF : " << neibhourLeft[i] ;
    //     else
    //         cout << " *** Corner left *** ";
    //     cout <<endl;

    //     cout << " Right : " ;
    //     if(neibhourRight[i])
    //         cout << x_pos[neibhourRight[i]] << " , " << y_pos[neibhourRight[i]] << " DOF : " << neibhourRight[i];
    //     else
    //         cout << " *** Corner Right *** ";
    //     cout <<endl;

    //     cout << " Top : " ;
    //     if(neibhourTop[i])
    //         cout << x_pos[neibhourTop[i]] << " , " << y_pos[neibhourTop[i]] << " DOF : " << neibhourTop[i];
    //     else
    //         cout << " *** Corner Top *** ";
    //     cout <<endl;

    //     cout << " Bottom : " ;
    //     if(neibhourBottom[i])
    //         cout << x_pos[neibhourBottom[i]] << " , " << y_pos[neibhourBottom[i]]<< " DOF : " << neibhourBottom[i] ;
    //     else
    //         cout << " *** Corner Bottom *** ";
    //     cout <<endl;
    // }


}

double *FTLE::obtainVelocity(double* uVelocityatPoint, double* vVelocityatPoint)
{
    double *uVector = VelocityVectFunction->GetComponent(0)->GetValues();
    double *vVector = VelocityVectFunction->GetComponent(1)->GetValues();

    for (int i = 0; i < N_Particles; i++)
    {
        if(!particleInsideDomain[i]) continue; // If particle is not inside the domain , then proceed to next particle

        int currentCellNo = cellIds[i];
        std::vector<int> neibhours = DOF_Neibhours_depth[currentCellNo];
        bool insideDomain = false;
        TBaseCell* cell;
        for ( int cellId = 0 ; cellId < neibhours.size(); cellId++)
        {
            cell =  fespace->GetCollection()->GetCell(cellId);

            if(cell->PointInCell(x_pos[i],y_pos[i]))
            {
                insideDomain = true;
                currentCellNo =  cellId;
                break;
            }
        }

        if(!insideDomain)  
        {
            particleInsideDomain[i] = false;
            continue;
        }
        double uVal[3];
        double vVal[3];
       VelocityVectFunction->GetComponent(0)->FindGradientLocal(cell,currentCellNo ,x_pos[i],y_pos[i],uVal);
       VelocityVectFunction->GetComponent(1)->FindGradientLocal(cell,currentCellNo ,x_pos[i],y_pos[i],vVal);
        
        // The first value will be the velocity value, others are derivative values
        uVelocityatPoint[i] = uVal[0];
        vVelocityatPoint[i] = vVal[0];

    }
}