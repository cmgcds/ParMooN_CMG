
/** ************************************************************************ 
* @brief     source file for ALE implementation in 3D
* @author    Thivin Anandh ( thivinanandh@iisc.ac.in )
* @date       02-Apr-2020
* @history    Final Working version for the Hehaheadral Cells  ( Conforming )
 ************************************************************************  */

#include<FE3D_ALE.h>
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <BoundFace.h>
#include <BoundComp3D.h>
#include <Assemble3D.h>
#include <AssembleMat3D.h>
#include <DirectSolver.h>
#include <FE3D.h>
#include <FESpace3D.h>
#include <Output3D.h>
#include <DirectSolver.h>
#include <string.h>
#include <math.h>
#include <ctime>
#include <LinAlg.h>
#include <FEVectFunct3D.h>
#include <sys/stat.h>
#include <vector>
#include <exception>
#include <fstream>
#include <algorithm>
#include <LinAlg.h>
#include <list>
#include <map>
#include <tuple>
#include <FEDatabase3D.h>
#include <BoundFace.h>
#include <IsoJointEqN.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <FEFunction3D.h>
#include <InterfaceJoint3D.h>
#include <NodalFunctional3D.h>
#include <mkl.h>
#include <bits/stdc++.h>
#ifdef  _CUDA
	#include<cuda.h>
	#include "helper_cuda.h"
	#include "nvToolsExt.h"
#endif  //_CUDA

#include <omp.h>
#include <chrono> 


const uint32_t colors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000, 0xffffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

#define PUSH_RANGE(name,cid) { \
    int color_id = cid; \
    color_id = color_id%num_colors;\
    nvtxEventAttributes_t eventAttrib = {0}; \
    eventAttrib.version = NVTX_VERSION; \
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
    eventAttrib.colorType = NVTX_COLOR_ARGB; \
    eventAttrib.color = colors[color_id]; \
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
    eventAttrib.message.ascii = name; \
    nvtxRangePushEx(&eventAttrib); \
}
#define POP_RANGE nvtxRangePop();

using namespace std::chrono; 
// 

// CONSTRUCTOR 
FE3D_ALE::FE3D_ALE(TCollection* _coll,	BoundCondFunct3D **_GridBoundaryConditions, BoundValueFunct3D **_GridBoundValues, 
					TFESpace3D* Meshfespace,TFEVectFunct3D* MeshVelocityVectFunction3D )
{
    GridBoundaryConditions  = _GridBoundaryConditions;
    GridBoundValues  = _GridBoundValues;
    coll              = _coll;

    //Get and Store all the Boundary Edges with  FREESLIP Boundary Condition

	sqstructure = new TSquareStructure3D(Meshfespace);
	sqstructure -> Sort();
	TSquareMatrix3D *SQMATRICES_GRID[9];
	double *RHS[3];
	double *Entries[9];
	int *GridKCol, *GridRowPtr;

	SqmatrixG11 = new TSquareMatrix3D(sqstructure);
	SqmatrixG12 = new TSquareMatrix3D(sqstructure);
	SqmatrixG13 = new TSquareMatrix3D(sqstructure);
	SqmatrixG21 = new TSquareMatrix3D(sqstructure);
	SqmatrixG22 = new TSquareMatrix3D(sqstructure);
	SqmatrixG23 = new TSquareMatrix3D(sqstructure);
	SqmatrixG31 = new TSquareMatrix3D(sqstructure);
	SqmatrixG32 = new TSquareMatrix3D(sqstructure);
	SqmatrixG33 = new TSquareMatrix3D(sqstructure);


	// get the Minimum Diameter of the cell collection
	InitialMeshMinDiameter = coll->GetShortestEdgeinCollection();
	InitialMeshMinVolume    = coll->GetminVolumeCollection();

	// If the cuda Solver Parameter is Set INitialise the cudaSolver Class
	#ifdef  _CUDA
	if(TDatabase::ParamDB->CUDASOLVERFLAG)
	{
		if(0 == (strcmp(TDatabase::ParamDB->CUDASOLVERTYPE, "REFACTORLU")))
		{
			cudaSolverReFactLU = new cudaRefactor();
		}
		else if(0 == (strcmp(TDatabase::ParamDB->CUDASOLVERTYPE, "REFACTORQR")))
		{
			cudaSolverLowLevQR = new cudaLowLevelQR();
		}

		else if (0 == (strcmp(TDatabase::ParamDB->CUDASOLVERTYPE, "REFACTORQR_OPTIMIZED")))
		{
         	cudaSolverLowLevQR_opt =  new cudaLowLevelQR_Optimised();	
		}
		else
		{
			cout << " ERROR : UNknown CUDASOLVERTYPE Selected in ParamFile " <<endl;
			cout << " INFO : Class : FE3D_ALE , Function : COnstructor , FIle : src/FE/FE3D_ALE.c"<<endl;
			cout << " Exiting Program " << endl;
			exit(0);     // THIVIN - EXIT Statement
  		}
		
	}
	
	// Get the nodal cordinates value and send them to the Device (GPU)
	// These arrays will provide mapping to local ref co-ordinates (xi, eta, zeta) given the local node numbering
	//### WARNING THIVIN :  It is assumed that all the cells in the mesh have same "Reftrans mapping" , 
	//                       Same Shape and same FE Order

	TBaseCell* t_cell = coll->GetCell(0);
	FE3D m_h_FEId = Meshfespace->GetFE3D(0, t_cell);
	RefTrans3D m_h_refTrans =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(m_h_FEId);


	int VertexPerCell = 8;         /// THIVIN -- HARDCODED VALUE
	m_h_N_VertexPerCell = VertexPerCell;


	TFE3D* FE_Obj_mesh = TFEDatabase3D::GetFE3D(m_h_FEId);
	TBaseFunct3D* bf_mesh = FE_Obj_mesh->GetBaseFunct3D();

	int N_BaseFunct_mesh = bf_mesh->GetDimension();
	int N_Points;

	//Crteate Cuda Stream 
	checkCudaErrors(cudaStreamCreate(&FE3d_stream));
	checkCudaErrors(cudaStreamCreate(&FE3d_stream_slip));

	checkCudaErrors(cudaStreamCreateWithFlags(&FE3d_stream,cudaStreamNonBlocking));
	checkCudaErrors(cudaStreamCreateWithFlags(&FE3d_stream_slip,cudaStreamNonBlocking));




	checkCudaErrors(cudaMalloc((void**)&m_d_N_VertexperCell,sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&m_d_N_NodalRefValues, sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&m_d_nodalFunctionalRefValues,3*N_BaseFunct_mesh*sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&m_d_nodalFunctionalRefValues_slip,3*N_BaseFunct_mesh*sizeof(double)));



	double *t_h_nodalValues_1, *t_h_nodalValues_2, *t_h_nodalValues_3;

	TNodalFunctional3D *nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(m_h_FEId);
	nf->GetPointsForAll(N_Points, t_h_nodalValues_1, t_h_nodalValues_2, t_h_nodalValues_3);

	checkCudaErrors(cudaMemcpyAsync(m_d_nodalFunctionalRefValues					 , t_h_nodalValues_1,sizeof(double)*N_BaseFunct_mesh,cudaMemcpyHostToDevice,FE3d_stream));
	checkCudaErrors(cudaMemcpyAsync(m_d_nodalFunctionalRefValues + N_BaseFunct_mesh   , t_h_nodalValues_2,sizeof(double)*N_BaseFunct_mesh,cudaMemcpyHostToDevice,FE3d_stream));
	checkCudaErrors(cudaMemcpyAsync(m_d_nodalFunctionalRefValues + 2*N_BaseFunct_mesh , t_h_nodalValues_3,sizeof(double)*N_BaseFunct_mesh,cudaMemcpyHostToDevice,FE3d_stream));


	checkCudaErrors(cudaMemcpyAsync(m_d_nodalFunctionalRefValues_slip					 , t_h_nodalValues_1,sizeof(double)*N_BaseFunct_mesh,cudaMemcpyHostToDevice,FE3d_stream_slip));
	checkCudaErrors(cudaMemcpyAsync(m_d_nodalFunctionalRefValues_slip + N_BaseFunct_mesh   , t_h_nodalValues_2,sizeof(double)*N_BaseFunct_mesh,cudaMemcpyHostToDevice,FE3d_stream_slip));
	checkCudaErrors(cudaMemcpyAsync(m_d_nodalFunctionalRefValues_slip + 2*N_BaseFunct_mesh , t_h_nodalValues_3,sizeof(double)*N_BaseFunct_mesh,cudaMemcpyHostToDevice,FE3d_stream_slip));


	#endif  //_CUDA


}


// TFEVectFunct3Dp
// ---------------------------- DOF PICKING ------------------------------------- //
void FE3D_ALE::pickDOFsOfFreeSlipBoundaries(TFESpace3D* gridfespace, std::vector<int> freeSlipBoundIds,std::vector<int> boundIds)
{
    int N_cells = coll->GetN_Cells();
    TBaseCell* currentCell;	
    int MaxLen;
    int N_Joints;
    const int* TmpLen; const int* TmpFV;
    BoundCond Bdcond;
    TBoundFace* Bdface;					// Pointer to Boundary Face in 3D Cell
	TBoundComp *BoundComp;
    TVertex* currentVertex;
    bool cell_setflag = false;
    int* GlobalNumbers;
    int* BeginIndex;
    int N_Movfaces = 0;
    GlobalNumbers = gridfespace->GetGlobalNumbers();
	BeginIndex = gridfespace->GetBeginIndex();	

	cout << " Function - Pick Free SLIP BD's " <<endl;
    // Member Variables  -- Do  not Initalise them 
    N_bd_FreeSlip_Vertex = 0;
    N_bd_FreeSlip_Cells = 0;
    N_bd_FreeSlip_Joints = 0;
	N_bd_EdgeFreeSlip_Vertex = 0;
    N_bd_EdgeFreeSlip_Cells = 0;
    N_bd_EdgeFreeSlip_Joints = 0;

	std::map<int,std::vector<int>> BdIdforDOF;
	//collect the number of Cells in Boundary 

	int N_ActiveBounds =  gridfespace->GetActiveBound();

	std::vector<int> dirchletDOFs;


	for( int cellId = 0 ; cellId < N_cells ; cellId++)
	{
       	currentCell = coll->GetCell(cellId);
		int* GlobalDOF = GlobalNumbers + BeginIndex[cellId];
		FE3D elementId = gridfespace->GetFE3D(cellId, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        N_Joints = currentCell->GetN_Joints();

		for ( int jointId = 0 ; jointId < N_Joints ; jointId++) 
		{
			TJoint* Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				Bdface = (TBoundFace*)Joint;
				BoundComp = Bdface->GetBoundComp();
	  			int bdid = BoundComp->GetID();
				if(bdid == 0)
				{
					int N_Vertices = TmpLen[jointId];	
					int *JointDOF = fedesc->GetJointDOF(jointId);				
					
					for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
					{
						int local_vertex_no = TmpFV[jointId*MaxLen + vert];
						currentVertex 	= currentCell->GetVertex(vert) ;
						int glob_vertex_no = GlobalDOF[JointDOF[vert]];
						dirchletDOFs.push_back(glob_vertex_no);
					}
			
				}
			}
		}
					
	}



	for (int cellNr = 0 ; cellNr < N_cells ; cellNr++)
    {
        currentCell = coll->GetCell(cellNr);
		int* GlobalDOF = GlobalNumbers + BeginIndex[cellNr];
		FE3D elementId = gridfespace->GetFE3D(cellNr, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        N_Joints = currentCell->GetN_Joints();

		for ( int jointId = 0 ; jointId < N_Joints ; jointId++)  // Joints refer to Faces in 3D
		{
			TJoint* Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				Bdface = (TBoundFace*)Joint;
				BoundComp = Bdface->GetBoundComp();
				int bdid = BoundComp->GetID();
				std::vector<int>::iterator it = std::find(freeSlipBoundIds.begin(), freeSlipBoundIds.end(), bdid);
				if(it != freeSlipBoundIds.end())                        
				{
					N_freeSurfaceJoints++;
					// if(cell_setflag == FALSE){
						N_freeSurfaceCells++;					
					// }   
					int *JointDOF = fedesc->GetJointDOF(jointId);
					int N_Vertices = TmpLen[jointId];
					N_bd_FreeSlip_Joints++;

					//freeSurfaceJoints.emplace_back(jointId);
					for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
					{
						int local_vertex_no = JointDOF[vert] ; // TmpFV[jointId*MaxLen + vert];
						currentVertex = currentCell->GetVertex(vert) ;
						int glob_vertex_no = GlobalDOF[JointDOF[vert]];
						if(N_bd_FreeSlip_Vertex == 0){
							std::vector<int> a;
							a.push_back(bdid*1000);
							a.push_back(local_vertex_no);
							a.push_back(cellNr);
							a.push_back(jointId);
							BdIdforDOF.insert(std::make_pair(glob_vertex_no,a));
							N_bd_FreeSlip_Vertex++;
						}
						// if (std::find(freeSurfaceVertex.begin(), freeSurfaceVertex.end(), glob_vertex_no) == freeSurfaceVertex.end()) 
						// {
							auto it = BdIdforDOF.find(glob_vertex_no);
							if ( it == BdIdforDOF.end() )   // First Time Entry of the Element
							{	
								std::vector<int> a;
								a.push_back(bdid*1000);
								a.push_back(local_vertex_no);
								a.push_back(cellNr);
								a.push_back(jointId);
								BdIdforDOF.insert(std::make_pair(glob_vertex_no,a));
								N_bd_FreeSlip_Vertex++;
							} 
							
							else if( it->second[0] != 1000*bdid )
							{
								it->second[0] = 9999;
								it->second.push_back(jointId);
								N_bd_FreeSlip_Vertex--;
							}
						// }

					}		
			
				}

			}

		}
	}
	


	N_bd_EdgeFreeSlip_Vertex = 0;
	N_bd_FreeSlip_Vertex = 0;
	for (auto const& pair: BdIdforDOF)
	{
		// auto it = std::find(dirchletDOFs.begin(),dirchletDOFs.end(),pair.first) ;
		if(pair.first < N_ActiveBounds)
		{
			if(pair.second[0] == 9999)
			{
				Bd_EdgeFreeSlip_Vertex.push_back(pair.first);
				Bd_EdgeFreeSlip_VertexLocal.emplace_back(pair.second[1]);
				Bd_EdgeFreeSlip_Cells.emplace_back(pair.second[2]);
				Bd_EdgeFreeSlip_Joints1.emplace_back(pair.second[3]);
				Bd_EdgeFreeSlip_Joints2.emplace_back(pair.second[4]);
				N_bd_EdgeFreeSlip_Vertex++;
			}
			else
			{
				Bd_FreeSlip_Vertex.push_back(pair.first);
				Bd_FreeSlip_VertexLocal.emplace_back(pair.second[1]);
				Bd_FreeSlip_Cells.emplace_back(pair.second[2]);
				Bd_FreeSlip_Joints.emplace_back(pair.second[3]);
				N_bd_FreeSlip_Vertex++;
			}
		}
    }



	cout << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  " <<endl;
		cout << " N FREE SLIP  : "<<N_bd_FreeSlip_Vertex<<endl;

	// for ( int i = 0 ; i < N_bd_FreeSlip_Vertex ; i++){
	// 	double x, y , z;
	// 	fespace_alemeshMovement->GetDOFPosition(Bd_FreeSlip_Vertex[i],x,y,z);
	// 	cout << Bd_FreeSlip_Vertex[i] << "(" << Bd_FreeSlip_Joints[i] <<")     " <<  " co ord : (" <<x<<", "<<y<<", "<<z<<" )         ;;;   "  ;
	// }
	// cout<<endl;

	// for ( int i = 0 ; i < N_bd_EdgeFreeSlip_Vertex ; i++)
	// 	cout << Bd_EdgeFreeSlip_Vertex[i] << "(" << Bd_EdgeFreeSlip_Joints1[i] <<","<< Bd_EdgeFreeSlip_Joints2[i] <<")    ";
	// cout<<endl;

	cout << " N FREESLIP EDGE : "<<N_bd_EdgeFreeSlip_Vertex<<endl;
	// cout << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  " <<endl;


	//Clear the Map Data Structiure
	BdIdforDOF.clear();


	// Allocate memory to the NOrmal Vectors
		// Resize the normal Vectors based on the number of free surface vertex available
	Bd_edge_normalA_1.resize(N_bd_EdgeFreeSlip_Vertex,0);
	Bd_edge_normalA_2.resize(N_bd_EdgeFreeSlip_Vertex,0);
	Bd_edge_normalA_3.resize(N_bd_EdgeFreeSlip_Vertex,0);

	Bd_edge_normalB_1.resize(N_bd_EdgeFreeSlip_Vertex,0);
	Bd_edge_normalB_2.resize(N_bd_EdgeFreeSlip_Vertex,0);
	Bd_edge_normalB_3.resize(N_bd_EdgeFreeSlip_Vertex,0);

	Bd_edge_TangentA_1.resize(N_bd_EdgeFreeSlip_Vertex,0);
	Bd_edge_TangentA_2.resize(N_bd_EdgeFreeSlip_Vertex,0);
	Bd_edge_TangentA_3.resize(N_bd_EdgeFreeSlip_Vertex,0);

	// Resize the normal Vectors based on the number of free surface vertex available
	Bd_normal_1.resize(N_bd_FreeSlip_Vertex,0);
	Bd_normal_2.resize(N_bd_FreeSlip_Vertex,0);
	Bd_normal_3.resize(N_bd_FreeSlip_Vertex,0);

	Bd_TangentA_1.resize(N_bd_FreeSlip_Vertex,0);
	Bd_TangentA_2.resize(N_bd_FreeSlip_Vertex,0);
	Bd_TangentA_3.resize(N_bd_FreeSlip_Vertex,0);

	Bd_TangentB_1.resize(N_bd_FreeSlip_Vertex,0);
	Bd_TangentB_2.resize(N_bd_FreeSlip_Vertex,0);
	Bd_TangentB_3.resize(N_bd_FreeSlip_Vertex,0);


	#ifdef  _CUDA
	//Memory Allocation Slip
	checkCudaErrors(cudaMalloc((void**)&m_d_Free_Slip_CellMapping,sizeof(int)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_Free_Slip_Vertex,sizeof(int)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_Free_Slip_Cells,sizeof(int)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_free_Slip_Joints,sizeof(int)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_free_Slip_VertexLocal,sizeof(int)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_free_Slip_Normal_1,sizeof(double)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_free_Slip_Normal_2,sizeof(double)*N_bd_FreeSlip_Vertex));
   checkCudaErrors(cudaMalloc((void**)&m_d_free_Slip_Normal_3,sizeof(double)*N_bd_FreeSlip_Vertex));  




	size_t N = sizeof(int)*N_bd_FreeSlip_Vertex;

	checkCudaErrors(cudaMemcpyAsync(m_d_Free_Slip_Vertex      ,Bd_FreeSlip_Vertex.data()      ,N,cudaMemcpyHostToDevice,FE3d_stream_slip));
   	checkCudaErrors(cudaMemcpyAsync(m_d_Free_Slip_Cells       ,Bd_FreeSlip_Cells.data()       ,N,cudaMemcpyHostToDevice,FE3d_stream_slip));
   	checkCudaErrors(cudaMemcpyAsync(m_d_free_Slip_VertexLocal ,Bd_FreeSlip_VertexLocal.data() ,N,cudaMemcpyHostToDevice,FE3d_stream_slip));
	checkCudaErrors(cudaMemcpyAsync(m_d_free_Slip_Joints      ,Bd_FreeSlip_Joints.data()      ,N,cudaMemcpyHostToDevice,FE3d_stream_slip));



	uniqueCellNumbers_freeSlip.resize(N_bd_FreeSlip_Vertex);
	m_h_Free_Slip_CellMapping.resize(N_bd_FreeSlip_Vertex);
	memcpy(uniqueCellNumbers_freeSlip.data(),Bd_FreeSlip_Cells.data(),N_bd_FreeSlip_Vertex*sizeof(int));

	std::sort(uniqueCellNumbers_freeSlip.begin(),uniqueCellNumbers_freeSlip.end());

	std::vector<int>::iterator iter;

   	iter = std::unique(uniqueCellNumbers_freeSlip.begin(),uniqueCellNumbers_freeSlip.end());

   	uniqueCellNumbers_freeSlip.resize(std::distance(uniqueCellNumbers_freeSlip.begin(),iter));


	for (int i = 0 ; i < N_bd_FreeSlip_Vertex ; i++)
   {
      iter = std::find(uniqueCellNumbers_freeSlip.begin(),uniqueCellNumbers_freeSlip.end(),Bd_FreeSlip_Cells[i]);
      if(iter != uniqueCellNumbers_freeSlip.end())
      {
		m_h_Free_Slip_CellMapping[i] =  iter - uniqueCellNumbers_freeSlip.begin();
      }
      
      else
      {
		cout << " ERROR : Could not find Cell NUmber " << Bd_FreeSlip_Cells[i]  << " i val : " << i <<endl;
		cout << " INFO  :  Function -  pick_freeSurface_SLIP_DOFs() , Class : FE3D_ALE  " <<endl;
		exit(0);
      }
      
   }

   checkCudaErrors(cudaMemcpyAsync(m_d_Free_Slip_CellMapping,m_h_Free_Slip_CellMapping.data(),sizeof(int)*N_bd_FreeSlip_Vertex,cudaMemcpyHostToDevice,FE3d_stream_slip));

	checkCudaErrors(cudaMalloc((void**)&X_cord_slip,sizeof(double)*8*uniqueCellNumbers_freeSlip.size()));
	checkCudaErrors(cudaMalloc((void**)&Y_cord_slip,sizeof(double)*8*uniqueCellNumbers_freeSlip.size()));
	checkCudaErrors(cudaMalloc((void**)&Z_cord_slip,sizeof(double)*8*uniqueCellNumbers_freeSlip.size()));

	
	
	#endif  //_CUDA

	

}	

//
void FE3D_ALE::Pick_free_surface_DOFs(TFESpace3D* _fespace, std::vector<int> BoundIds,TCollection* coll)
{
	N_freeSurfaceJoints = 0;
	N_freeSurfaceCells = 0;
	N_freeSurfaceVertex = 0 ;
	fespace_alemeshMovement = _fespace;
	
	// -------------------- Variable Declaraions -------------------------//
		
    int* N_Movfaces; 					// Local Face id's of the boundary Domain 
	TBaseCell* currentCell;   			// Cell pointer for the current cell
	TJoint* Joint;						// Pointer to Joint in 3D Cell
	TBoundFace* Bdface;					// Pointer to Boundary Face in 3D Cell
	TBoundComp *BoundComp;
	TVertex* currentVertex;


	int N_Cells = coll->GetN_Cells();   // N_cells in given FE system
	int N_Vertices; 
	int N_Joints;
	int TMP_DOF = 0;
	const int *TmpFV,*TmpLen;
	int MaxLen;
	int *GlobalDOF;
	int *BeginIndex;
	int* GlobalNumbers;

	GlobalNumbers = fespace_alemeshMovement->GetGlobalNumbers();
	BeginIndex = fespace_alemeshMovement->GetBeginIndex();	

	// for ( int k = 0 ; k < N_Cells*8 ; k++)
	// {
	// 	cout << "loc : " << k <<  "  Glob  " <<   GlobalNumbers[k] << "  " <<endl;
 	// }
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////// IDEntify the Edge DOF's in the given Free Surface and Exclude them from selection in Free Surface DOF ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Identify the Edge DOF's in the Free Surface

	std::map<int,int> GlobalBdIdforDOF;

	for( int cellId = 0 ; cellId < N_Cells ; cellId++)
	{
		currentCell = coll->GetCell(cellId); 	    // Obatin the pointer to current cell from Collection
		currentCell->SetGlobalCellNo(cellId);
		GlobalDOF = GlobalNumbers + BeginIndex[cellId];
								
		FE3D elementId = fespace_alemeshMovement->GetFE3D(cellId, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        N_Joints = currentCell->GetN_Joints();

		for ( int jointId = 0 ; jointId < N_Joints ; jointId++) 
		{
			Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				Bdface = (TBoundFace*)Joint;
				BoundComp = Bdface->GetBoundComp();
				int bdid = BoundComp->GetID();
				int N_Vertices = TmpLen[jointId];	
				int *JointDOF = fedesc->GetJointDOF(jointId);				
			
				for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
				{
					int local_vertex_no =   JointDOF[vert]   ; //  TmpFV[jointId*MaxLen + vert];
					currentVertex = currentCell->GetVertex(vert) ;
					int glob_vertex_no = GlobalDOF[JointDOF[vert]];
					if (GlobalBdIdforDOF[glob_vertex_no] != 1000 + bdid + 1 )
						GlobalBdIdforDOF[glob_vertex_no] += 1000 + bdid + 1;
				}
			}
		}
					
	}

	std::vector<int> GlobaledgeDOFs;
	for (auto it : GlobalBdIdforDOF) 
	{
		if (it.second > 2000)  // Edge Nodes 
		{
			GlobaledgeDOFs.push_back(it.first);
		}
	}

	GlobalBdIdforDOF.clear();

	for( int cellId = 0 ; cellId < N_Cells ; cellId++)
	{
		currentCell = coll->GetCell(cellId); 	    // Obatin the pointer to current cell from Collection
		currentCell->SetGlobalCellNo(cellId);
		GlobalDOF = GlobalNumbers + BeginIndex[cellId];
								
		FE3D elementId = fespace_alemeshMovement->GetFE3D(cellId, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        N_Joints = currentCell->GetN_Joints();

		for ( int jointId = 0 ; jointId < N_Joints ; jointId++) 
		{
			Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				Bdface = (TBoundFace*)Joint;
				BoundComp = Bdface->GetBoundComp();
				int bdid = BoundComp->GetID();
				
				if (std::find(BoundIds.begin(), BoundIds.end(), bdid) != BoundIds.end())
				{
					int N_Vertices = TmpLen[jointId];	
					int *JointDOF = fedesc->GetJointDOF(jointId);				
				
					for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
					{
						int local_vertex_no =  JointDOF[vert] ; //TmpFV[jointId*MaxLen + vert];
						currentVertex = currentCell->GetVertex(vert) ;
						int glob_vertex_no = GlobalDOF[JointDOF[vert]];
						GlobalBdIdforDOF[glob_vertex_no] += 1000 + bdid + 1;
					}
				}
			}
		}
					
	}



	std::vector<int> freeSurfDOFs;
	for (auto it : GlobalBdIdforDOF) 
	{
		freeSurfDOFs.push_back(it.first);
	}

	std::sort(GlobaledgeDOFs.begin(),GlobaledgeDOFs.end());
	std::sort(freeSurfDOFs.begin(),freeSurfDOFs.end());

	FreeSurfaceEdgeDOFs.resize(freeSurfDOFs.size()); 
	std::vector<int>::iterator it, st;
	it = set_intersection(GlobaledgeDOFs.begin(), 
                          GlobaledgeDOFs.end(), 
                          freeSurfDOFs.begin(), 
                          freeSurfDOFs.end(), 
                          FreeSurfaceEdgeDOFs.begin()); 

	FreeSurfaceEdgeDOFs.resize(it -FreeSurfaceEdgeDOFs.begin() );
	cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" <<endl;
	cout << " size of fresurf edge: " << FreeSurfaceEdgeDOFs.size()<<endl;
	cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" <<endl;


	////////////////////////////// --------------- END------------------- ///////////////////////////////////////////////////////////////////////////////

	// MovBoundVert = new TVertex*[300];
	int counter = 0;
	// Loop over all the cells to collect the DATA VALUES
	for( int cellId = 0 ; cellId < N_Cells ; cellId++)
	{
		currentCell = coll->GetCell(cellId); 	    // Obatin the pointer to current cell from Collection
		currentCell->SetGlobalCellNo(cellId);
		GlobalDOF = GlobalNumbers + BeginIndex[cellId];
							
		FE3D elementId = fespace_alemeshMovement->GetFE3D(cellId, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
	
		N_Joints = currentCell->GetN_Joints();
		bool cell_setflag = FALSE;

		currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

		// Obtain the Vertex and face information of the current cell
		// Loop Over the Joints in the cell 
		for ( int jointId = 0 ; jointId < N_Joints ; jointId++)  // Joints refer to Faces in 3D
		{
			Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				// cout << "  Boundary Face " <<endl;
				Bdface = (TBoundFace*)Joint;
	  			BoundComp = Bdface->GetBoundComp();
	  			int bdid = BoundComp->GetID();
				std::vector<int>::iterator it = std::find(BoundIds.begin(), BoundIds.end(), bdid);
				if(it != BoundIds.end())                     // CHANGE THIS HARDCODED VALUE
				{
					N_freeSurfaceJoints++;
					if(cell_setflag == FALSE){
						N_freeSurfaceCells++;					
						cell_setflag = TRUE ;
					}   
					int *JointDOF = fedesc->GetJointDOF(jointId);
					N_Vertices = TmpLen[jointId];
					N_Movfaces++;
					//freeSurfaceJoints.emplace_back(jointId);
					for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
					{
						int local_vertex_no =  JointDOF[vert]  ;   // TmpFV[jointId*MaxLen + vert];
						currentVertex = currentCell->GetVertex(vert) ;
						int glob_vertex_no = GlobalDOF[JointDOF[vert]];
						if (std::find(freeSurfaceVertex.begin(), freeSurfaceVertex.end(), glob_vertex_no) == freeSurfaceVertex.end()) 
						{
							// if (std::find(edgeDOFs.begin(), edgeDOFs.end(), glob_vertex_no) == edgeDOFs.end()) 
							// {
								// MovBoundVert[N_freeSurfaceVertex] = currentVertex;
								freeSurfaceVertex.emplace_back(glob_vertex_no);
								freeSurfaceVertexLocal.emplace_back(local_vertex_no);
								freeSurfaceCells.emplace_back(cellId);
								freeSurfaceJoints.emplace_back(jointId);
								N_freeSurfaceVertex++;
							// }	
						}
						if(currentVertex->GetClipBoard() != -5 || N_freeSurfaceVertex == 0)
						{						
							currentVertex->SetClipBoard(-5);
							
						}
					}
				}
			}
		}
	}


	
	// cout << " ----- GLOBAL DOF's After remove " <<endl;
	// for ( int k = 0 ; k < freeSurfaceVertex.size() ; k++)
	// 	cout << " " << freeSurfaceVertex[k]  << "  Joint :  "<<freeSurfaceJoints[k] << "  Cells : "  << freeSurfaceCells[k] << endl;
	
	// cout << " ----------------------------" <<endl;

	// Collect all the Cells and Save them in Tvertex** POinter
	
	MovCells = new TBaseCell*[N_freeSurfaceVertex];
	MovJoints = new TJoint*[N_freeSurfaceJoints];


	for( int i = 0  ; i < N_freeSurfaceVertex ; i++){
		MovCells[i] = coll->GetCell(freeSurfaceCells[i]);
		// cout << freeSurfaceVertex[i] << "(" << freeSurfaceJoints[i] <<")     " ;
	}
	// cout <<endl;

	cout << " No of Free Surface DOF after DOF picking: " << N_freeSurfaceVertex<<endl;

	// ---------------------------------------------------------------------------------------------------------- // 
	// ------------ Allocation of memory for Free Srface DOF's Normal vectors and theor local Storages -----------//
	// ---------------------------------------------------------------------------------------------------------- // 

	xi_1.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face
	xi_2.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face 

	xi_freeSurfaceNodes.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face
	eta_freeSurfaceNodes.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face
	zeta_freeSurfaceNodes.resize(N_freeSurfaceVertex,0); // Arrays to store temp XI valu on 2d Face

	// Resize the normal Vectors based on the number of free surface vertex available
	meshVelocityatFreeSurface_1.resize(N_freeSurfaceVertex,0);
	meshVelocityatFreeSurface_2.resize(N_freeSurfaceVertex,0);
	meshVelocityatFreeSurface_3.resize(N_freeSurfaceVertex,0);

	Bd_FreeSurf_normal_1.resize(N_freeSurfaceVertex,0);
	Bd_FreeSurf_normal_2.resize(N_freeSurfaceVertex,0);
	Bd_FreeSurf_normal_3.resize(N_freeSurfaceVertex,0);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////// Allocate memory in Device for Normal Calculations            /////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	#ifdef  _CUDA	
	
	// This part is done Asynchronously <S o that the remaining Sequential part can be executed while the data transfer 
	// Takes place 

	checkCudaErrors(cudaMalloc((void**)&m_d_FreeSurfCellMapping,sizeof(int)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_FreeSurfaceVertex,sizeof(int)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_FreeSurfaceCells,sizeof(int)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_freeSurfaceJoints,sizeof(int)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_freeSurfaceVertexLocal,sizeof(int)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_freeSurfaceNormal_1,sizeof(double)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_freeSurfaceNormal_2,sizeof(double)*N_freeSurfaceVertex));
	checkCudaErrors(cudaMalloc((void**)&m_d_freeSurfaceNormal_3,sizeof(double)*N_freeSurfaceVertex));	

	// Allocate momory to the Device for X,Y,Z CoOrdinates
	double* m_d_Vertex_Xcoord;
	double* m_d_Vertex_Ycoord;
	double* m_d_Vertex_Zcoord;
	
	checkCudaErrors(cudaMalloc((void**)&m_d_Vertex_Zcoord,sizeof(double)*N_freeSurfaceVertex));	


	checkCudaErrors(cudaMalloc((void**)&m_d_N_FreeSurfaceVertex,sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&m_d_N_NodalRefValues,sizeof(int)));

	
	

	size_t N = sizeof(int)*N_freeSurfaceVertex;
	
	checkCudaErrors(cudaMemcpyAsync(m_d_FreeSurfaceVertex      ,freeSurfaceVertex.data()      ,N,cudaMemcpyHostToDevice,FE3d_stream));
	checkCudaErrors(cudaMemcpyAsync(m_d_FreeSurfaceCells       ,freeSurfaceCells.data()       ,N,cudaMemcpyHostToDevice,FE3d_stream));
	checkCudaErrors(cudaMemcpyAsync(m_d_freeSurfaceVertexLocal ,freeSurfaceVertexLocal.data() ,N,cudaMemcpyHostToDevice,FE3d_stream));
	checkCudaErrors(cudaMemcpyAsync(m_d_freeSurfaceJoints      ,freeSurfaceJoints.data()      ,N,cudaMemcpyHostToDevice,FE3d_stream));


	#endif  //_CUDA


	// ---------------------------------------------------------------------------------------------------------- // 
	// FIND THE LOCAL CO ORDINATES of the GLOBAL DOF's and SAve them in an Array for Furture Normal Calculationn //
	// ---------------------------------------------------------------------------------------------------------- // 

	BF3DRefElements RefElement ;	
	RefTrans3D RefTrans , RefTransVelocity;TRefTrans3D *F_K ;
	for ( int i = 0 ; i < N_freeSurfaceVertex ; i ++)
	{
		int vertex_number = i;
		
		int cellNr = freeSurfaceCells[i];
		currentCell =  coll->GetCell(cellNr);

		// Parameters for the Mesh Velocity FESPACE	---- /////
		FE3D FEId = fespace_alemeshMovement->GetFE3D(cellNr, currentCell);
		TFE3D* ele = TFEDatabase3D::GetFE3D(FEId);
		RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
		RefTrans3D referenceTransformation =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);		

		int* GlobalNumbers_Mesh = fespace_alemeshMovement->GetGlobalNumbers();
		int* BeginIndex_Mesh = fespace_alemeshMovement->GetBeginIndex();
		int* Numbers_Mesh = GlobalNumbers_Mesh + BeginIndex_Mesh[cellNr];


		int JointNumber = freeSurfaceJoints[i];

		switch(referenceTransformation)     // Reftrans of Velocity
    	{
			case HexaTrilinear: // HEXATRILINEAR 
			{	
				// RefTrans = HexaTrilinear;
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((THexaTrilinear *)F_K)->SetCell(currentCell);

				int localNodeNUmber = freeSurfaceVertexLocal[vertex_number];

				((THexaTrilinear *)F_K)->GetRefvaluesfromLocalNodeNumber(FEId, localNodeNUmber, xi_freeSurfaceNodes[vertex_number]
															, eta_freeSurfaceNodes[vertex_number], zeta_freeSurfaceNodes[vertex_number]);

				((THexaTrilinear *)F_K)->GetRefValuesfromJointid(JointNumber,xi_freeSurfaceNodes[vertex_number]
											,eta_freeSurfaceNodes[vertex_number],zeta_freeSurfaceNodes[vertex_number]
											,xi_1[vertex_number],xi_2[vertex_number]);
				break;
			}

			case TetraAffin :  
			{
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((TTetraAffin *)F_K)->SetCell(currentCell);

				int localNodeNUmber = freeSurfaceVertexLocal[vertex_number];

				((TTetraAffin *)F_K)->GetRefvaluesfromLocalNodeNumber(FEId, localNodeNUmber, xi_freeSurfaceNodes[vertex_number]
															, eta_freeSurfaceNodes[vertex_number], zeta_freeSurfaceNodes[vertex_number]);
				

				// The below Step is redundant , as the s and t ( xi_1 and xi_2 ) values will not be used for calculation of 
				// normals at the Tetraheadral Cell, However to maintain code consistency with hexalinear, we will assign the values to any of the 
				// twi reference co-ordinates 
				xi_1[vertex_number] = xi_freeSurfaceNodes[vertex_number];
				xi_2[vertex_number] = eta_freeSurfaceNodes[vertex_number];
				break;
			}


			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class FE3D_ALE , Function get_surface_normals " <<endl;
				exit(0);
			}
		}
	}


	#ifdef _CUDA
	// ---------------------------------------------------------------------------------------------------------- // 
	// ------------------------------ Find the Cell Mapping to be used within Device Memory -------------------- //
	// ---------------------------------------------------------------------------------------------------------- // 

	// all the Cells' xyz coordinates will be saved in a Single array , So we need to create a mapping of , which cell's co-ordinates start\
	// from which position in the Global 1d Array in Device. 
		
	//Get the unique Collection of Cell Numbers 
	uniqueCellNumbers.resize(N_freeSurfaceVertex);
	m_h_FreeSurfCellMapping.resize(N_freeSurfaceVertex);
	memcpy(uniqueCellNumbers.data(),freeSurfaceCells.data(),N_freeSurfaceVertex*sizeof(int));

	std::sort(uniqueCellNumbers.begin(),uniqueCellNumbers.end());
	std::vector<int>::iterator iter;

	iter = std::unique(uniqueCellNumbers.begin(),uniqueCellNumbers.end());

	uniqueCellNumbers.resize(std::distance(uniqueCellNumbers.begin(),iter));

	for (int i = 0 ; i < N_freeSurfaceVertex ; i++)
	{
		iter = std::find(uniqueCellNumbers.begin(),uniqueCellNumbers.end(),freeSurfaceCells[i]);
		if(iter != uniqueCellNumbers.end())
		{
			m_h_FreeSurfCellMapping[i] =  iter - uniqueCellNumbers.begin();
		}
		
		else
		{
			cout << " ERROR : Could not find Cell NUmber " << freeSurfaceCells[i]  << " i val : " << i <<endl;
			cout << " INFO  :  Function -  pick_freeSurface_DOFs() , Class : FE3D_ALE  " <<endl;
			exit(0);
		}
		
	}

	checkCudaErrors(cudaMemcpyAsync(m_d_FreeSurfCellMapping,m_h_FreeSurfCellMapping.data(),sizeof(int)*N_freeSurfaceVertex,cudaMemcpyHostToDevice,FE3d_stream));
	// checkCudaErrors(cudaMalloc((void**)&m_d_Vertex_Xcoord,sizeof(double)*8*uniqueCellNumbers.size()));
	// checkCudaErrors(cudaMalloc((void**)&m_d_Vertex_Ycoord,sizeof(double)*8*uniqueCellNumbers.size()));
	// checkCudaErrors(cudaMalloc((void**)&m_d_Vertex_Zcoord,sizeof(double)*8*uniqueCellNumbers.size()));

	checkCudaErrors(cudaMalloc((void**)&X_cord,sizeof(double)*m_h_N_VertexPerCell*uniqueCellNumbers.size()));
	checkCudaErrors(cudaMalloc((void**)&Y_cord,sizeof(double)*m_h_N_VertexPerCell*uniqueCellNumbers.size()));
	checkCudaErrors(cudaMalloc((void**)&Z_cord,sizeof(double)*m_h_N_VertexPerCell*uniqueCellNumbers.size()));

	
	#endif // _CUDA

}


// ------------------------- MESH ASSEMBLY MODULE ----------------------------------------//
void FE3D_ALE::GridCoeffs(int n_points, double *x, double *y,double *z, double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}


void FE3D_ALE::Assembly_poisson_3D(double quad_wt, double *coeff, double *param, double hK, 
									double **derivatives, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double *N = derivatives[0], *Nx = derivatives[1], *Ny = derivatives[2] , *Nz = derivatives[3];
  double **K11, **K12, **K13, **K21, **K22, **K23 , **K31, **K32, **K33, *F1, *F2 , *F3;
  K11 = LocMatrices[0];
  K12 = LocMatrices[1];
  K13 = LocMatrices[2];
  K21 = LocMatrices[3];
  K22 = LocMatrices[4];
  K23 = LocMatrices[5];
  K31 = LocMatrices[6];
  K32 = LocMatrices[7];
  K33 = LocMatrices[8];
  
  F1 = LocRhs[0];
  F2 = LocRhs[1];
  F3 = LocRhs[1];

  for (int i = 0; i < N_BaseFuncts[0]; i++){
    for (int j = 0; j < N_BaseFuncts[0]; j++){
		K11[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j] + Nz[i]*Nz[j]);
		K22[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j] + Nz[i]*Nz[j]);
		K33[i][j] += quad_wt * (Nx[i]*Nx[j] + Ny[i]*Ny[j] + Nz[i]*Nz[j]);  

      //cout << "Nx[ " << i << "] = " << Nx[i] << endl;
      	K12[i][j] += 0.;
		K13[i][j] += 0.;
      	K21[j][i] += 0.;
		K23[j][i] += 0.;
		K31[j][i] += 0.;
		K32[j][i] += 0.;
		
      /* NON LINEAR PART */
    }
    /* RHS */

    F1[i] = 0.;
    F2[i] = 0.;
	F3[i] = 0.;
  }


}


void FE3D_ALE::Solve_mesh_velocity_interior(TFESpace3D* Meshfespace, TFEVectFunct3D* MeshVelocityVectFunction3D)
{
	char UString[] = "T";
	char NameString[] = "name";
	char CString[] = "C";

	int N_DOF = Meshfespace->GetN_DegreesOfFreedom();
	int N_Active = Meshfespace->GetActiveBound();
    
    // cout  << "N_Cells = " << coll -> GetN_Cells() << endl;
	// cout << "Degrees of Freedom = " << N_DOF  << "    N_Active = " << N_Active << endl;
	double *sol = new double[3*N_DOF]();
	double *rhs =  new double[3*N_DOF]();
	
	double* meshVelocityValues;
	meshVelocityValues = MeshVelocityVectFunction3D->GetValues();
    
    TAuxParam3D *Meshaux = new TAuxParam3D(1, 0, 0, 0, &Meshfespace, NULL, NULL, NULL, NULL, 0, NULL);
    
    int N_Terms = 4;    
  	int *SpacesNumbers = new int[N_Terms](); 
  	int N_Matrices = 9;    	
	int N_RHS = 3;   
	int *rowspace = new int[N_Matrices]();
	int *columnspace = new int[N_Matrices](); 
	int *rhsspace = new int[N_RHS]();
    
    
    MultiIndex3D AllDerivatives[4] = {D000, D100, D010,D001};
     
    TDiscreteForm3D* discreteform = new TDiscreteForm3D(UString, UString, N_Terms, AllDerivatives,
                                        SpacesNumbers, N_Matrices, N_RHS, rowspace, columnspace, rhsspace,
										                    Assembly_poisson_3D, GridCoeffs, NULL); 
    
    // --------------------- START OF MATRIX STRUCTURE DECLARATION -------------------//



	SQMATRICES_GRID[0] = SqmatrixG11;  SQMATRICES_GRID[1] = SqmatrixG12; SQMATRICES_GRID[2] = SqmatrixG13;
	SQMATRICES_GRID[3] = SqmatrixG21;  SQMATRICES_GRID[4] = SqmatrixG22; SQMATRICES_GRID[5] = SqmatrixG23;
	SQMATRICES_GRID[6] = SqmatrixG31;  SQMATRICES_GRID[7] = SqmatrixG32; SQMATRICES_GRID[8] = SqmatrixG33;
	

	RHS[0] = rhs; RHS[1] = rhs + N_DOF; RHS[2] = rhs + 2*N_DOF;

	Entries[0] = SqmatrixG11->GetEntries(); Entries[1] = SqmatrixG12->GetEntries(); Entries[2] = SqmatrixG13->GetEntries();
	Entries[3] = SqmatrixG21->GetEntries(); Entries[4] = SqmatrixG22->GetEntries(); Entries[5] = SqmatrixG23->GetEntries();
	Entries[6] = SqmatrixG31->GetEntries(); Entries[7] = SqmatrixG32->GetEntries(); Entries[8] = SqmatrixG33->GetEntries();
	
	GridKCol = sqstructure->GetKCol(); GridRowPtr = sqstructure->GetRowPtr();
	
	// ---------------- START OF ASSEMBLY 3D FUNCTION -----------//
	fesp[0] = Meshfespace;   // Type of FE Space to be used for Blocks in A Matrix
    ferhs[0] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
    ferhs[1] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	ferhs[2] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	
    TAssembleMat3D *MeshMatAssemble = new TAssembleMat3D(1, &Meshfespace, 9, SQMATRICES_GRID, 0, NULL,
														 3, RHS, ferhs, discreteform, GridBoundaryConditions,
														  GridBoundValues, Meshaux);
	MeshMatAssemble->Init();
	MeshMatAssemble->Reset();
	MeshMatAssemble->Assemble3D();

	// N_Active = Non Dirichlet DOF's
	int N_BDDof = N_DOF - N_Active;

	int N_ActiveBoundary = SqmatrixG11->GetActiveBound();
	

	const int N_DOF_perBlock =  SqmatrixG11->GetN_Rows();
	const int N_Active_perBlock = SqmatrixG11->GetActiveBound();
	const int N_Dirichlet_perBlock = N_DOF - N_Active;
	GridRowPtr = sqstructure->GetRowPtr();

	
	// // Memset the Antidiagonal arrays to be zero. 
	memset(SqmatrixG12->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] ) *sizeof(double));
	memset(SqmatrixG13->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG21->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG23->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG31->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG32->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));

	
	// Assign the u.n boundary Condition for the Free Slip Slides
	impose_FreeSlip_BoundaryCondition(sqstructure,SqmatrixG11->GetEntries(),SqmatrixG12->GetEntries(),SqmatrixG13->GetEntries(),SqmatrixG21->GetEntries(),SqmatrixG22->GetEntries()
									,SqmatrixG23->GetEntries(),SqmatrixG31->GetEntries(),SqmatrixG32->GetEntries(),SqmatrixG33->GetEntries(),rhs,N_DOF,N_ActiveBoundary);
	


	// Assign the Velocity Values to rhs from U
		// Transfer the Solution Calculated from "U" as a Dirichlet boundary to the Mesh Velocity Solution
	for ( int i = 0 ; i < N_freeSurfaceVertex ; i++)
	{	
		int globDOF 	= 	freeSurfaceVertex[i];
		std::vector<int>::iterator edgeDOF_it 	=   std::find(FreeSurfaceEdgeDOFs.begin(), FreeSurfaceEdgeDOFs.end(),globDOF ); 
		if ( edgeDOF_it != FreeSurfaceEdgeDOFs.end())  // Edge DOF
		{
			
			if(fabs(Bd_FreeSurf_normal_2[i]) > 1e3)
			{
				meshVelocityatFreeSurface_2[i] 	+=  meshVelocityatFreeSurface_1[i]*(Bd_FreeSurf_normal_1[i]/Bd_FreeSurf_normal_2[i]) +
													meshVelocityatFreeSurface_3[i]*(Bd_FreeSurf_normal_3[i]/Bd_FreeSurf_normal_2[i]) ;
			}

			rhs[globDOF           ] =  0;
			rhs[globDOF +  1*N_DOF] =  meshVelocityatFreeSurface_2[i];
			rhs[globDOF +  2*N_DOF] =  0;
		}

		else
		{
			if(fabs(Bd_FreeSurf_normal_2[i]) > 1e3)
			{
				meshVelocityatFreeSurface_2[i] 	+=  meshVelocityatFreeSurface_3[i]*(Bd_FreeSurf_normal_3[i]/Bd_FreeSurf_normal_2[i]) ;
			}

			rhs[globDOF           ] =  meshVelocityatFreeSurface_1[i];
			rhs[globDOF +  1*N_DOF] =  meshVelocityatFreeSurface_2[i];
			rhs[globDOF +  2*N_DOF] =  0;
		}
		

	}

	
	// ----------- DIRECT SOLVER -------------------------------- //
	cout << " RHSSS Velocity Solution Norm : " <<Ddot(3*N_DOF,rhs,rhs)<<endl;

	auto start = high_resolution_clock::now();
	PardisoDirectSolver_without_removing_dirichlet(SQMATRICES_GRID, 3, 3, sol, rhs);


	// Imposing Free Slip Boundary Condition
	// cout << "3*  NDOF " << 3*N_DOF <<endl;
 
	// auto end =  high_resolution_clock::now(); 
	// auto duration = duration_cast<milliseconds>(end - start);
	cout << " Mesh Velocity Solution Norm : " <<Ddot(3*N_DOF,sol,sol)<<endl;


	// -----------END---- DIRECT SOLVER -------------------------------- //
		// getchar();

	memcpy( meshVelocityValues, sol , SizeOfDouble*3*N_DOF);


	// Release the Matrix Storage Parameters
	for ( int i = 0 ; i < N_Matrices ; i++)
		SQMATRICES_GRID[i]->Reset();

    for (int i_rhs = 0; i_rhs < 3*N_DOF; i_rhs++)
      rhs[i_rhs] = 0;


	// std::cout << " Time Taken to Solve : "<<duration.count() << " ms"<<std::endl;
}


// ------------------------- GET SURFACE NORMALS ----------------------------------------//
void FE3D_ALE::get_surface_normals(TFEVectFunct3D* MeshVelo_FEvect,TFEVectFunct3D* Velocity_FEvect )
{
	//	variable Declarations
	// cout << " REACHED - Get surface Normals " <<endl;
	
	TBaseCell* currentCell;  TVertex* currentVertex;	
	FE3D FEId, FEId_velocity; TFE3D *ele , *ele_Velocity;
	BF3DRefElements RefElement , RefElement_Velocity;	
	RefTrans3D RefTrans , RefTransVelocity;TRefTrans3D *F_K , *F_K_V;
	TBaseFunct3D *bf;
	TFESpace3D *Velocity_FESpace,*Mesh_FESpace;
	TFEFunction3D*  Velocity_FEFunction[3];


	//  Get the Compunents of Velocity_FEvect and MeshVelo_fevect
	int meshVelComp = MeshVelo_FEvect->GetN_Components();
	int velComp = Velocity_FEvect->GetN_Components();

	//get FE Space for Velocity & Mesh
	Velocity_FESpace = Velocity_FEvect->GetFESpace3D();
	Mesh_FESpace = MeshVelo_FEvect->GetFESpace3D();	

	double* VelocityArray = Velocity_FEvect->GetComponent(0)->GetValues();
	double* MeshVelocityArray = MeshVelo_FEvect->GetComponent(0)->GetValues();

	// Create an FEFunction array for all the the FEVectfunction3d objects and to Store Value Arrays from the respective FEFunctions
	for ( int i = 0 ; i < velComp ; i++)
	{   
		Velocity_FEFunction[i] = Velocity_FEvect->GetComponent(i);
	}

	// Get the Lengths of the Value Arrays of Respective FEFunction Values
	int MeshVeloLength = MeshVelo_FEvect->GetComponent(0)->GetLength();
	int VelocityLength = Velocity_FEvect->GetComponent(0)->GetLength();

	// for ( int i = 0; i < VelocityLength ; i++)
	// {
	// 	VelocityArray[i                     ] = i;
	// 	VelocityArray[i + VelocityLength    ] = i + 100;
	// 	VelocityArray[i + (VelocityLength*2)] = i + 200;	
	// }

	double norm1 = 0.,norm2 = 0.,norm3 = 0.,normLen = 0.;

	//Calculate normal for Each Node in the Joint. 
	// Hard Code the Values of Joint id's for the given HEXAHEADRAL - TRILINEAR
	for ( int i = 0 ; i < N_freeSurfaceVertex ; i ++)  // N_freeSurfaceVertex
	{
		int vertex_number = i;
		// cout << " --------------------------- vertex number " <<i <<"---------------------- " <<endl;
		// cout << " Local DOF of the Vertex : " << freeSurfaceVertexLocal[i] ;
		//cout  << " Global DOF of the Vertex : " << freeSurfaceVertex[i] <<endl;
		int cellNr = freeSurfaceCells[i];
		currentCell =  coll->GetCell(cellNr);

		// Parameters for the Mesh Velocity FESPACE	---- /////
		FEId = Mesh_FESpace->GetFE3D(cellNr, currentCell);
		ele = TFEDatabase3D::GetFE3D(FEId);
		RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
		RefTrans3D referenceTransformation =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////


		// Parameters for the  Velocity FESPACE	---- /////
		FEId_velocity = Velocity_FESpace->GetFE3D(cellNr, currentCell);
		ele_Velocity = TFEDatabase3D::GetFE3D(FEId_velocity);
		RefElement_Velocity = TFEDatabase3D::GetRefElementFromFE3D(FEId_velocity);
		RefTrans3D RefTransVelocity =  ele_Velocity->GetRefTransID();
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////

		//----------- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //
		int* GlobalNumbers_Velocity = Velocity_FESpace->GetGlobalNumbers();
		int* BeginIndex_velocity = Velocity_FESpace->GetBeginIndex();
		int* Numbers_Velocity  = GlobalNumbers_Velocity + BeginIndex_velocity[cellNr];

		int* GlobalNumbers_Mesh = Mesh_FESpace->GetGlobalNumbers();
		int* BeginIndex_Mesh = Mesh_FESpace->GetBeginIndex();
		int* Numbers_Mesh = GlobalNumbers_Mesh + BeginIndex_Mesh[cellNr];
		//----END---- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //

		// SETUP the Values of s and T that needs to be sent to "GetOuterNormal"/ getTangentVectors to get the normal at the point
		double x_coord ,y_coord,z_coord;
		double xi_temp , eta_temp , zeta_temp;
		

		int JointNumber = freeSurfaceJoints[i];

		//////////// CODE BLOCK A2 - Calculate Joint Normal using "MESH" velocity FE Space  ////////////////////
		////// Note : This Code Block is Repeated Block of Code to calculate Normal of a Face at the mesh Face //////
		////// Note : This Code is for "MESH" Velocity FE Space only /////////////////
						
		switch(referenceTransformation)     // Reftrans of MESH Velocity
    	{
			case HexaTrilinear: // HEXATRILINEAR 
			{	
				//RefTrans = HexaTrilinear;
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((THexaTrilinear *)F_K)->SetCell(currentCell);
				( (THexaTrilinear *) F_K )->GetOuterNormal(JointNumber,xi_1[i],xi_2[i],
											Bd_FreeSurf_normal_1[i], Bd_FreeSurf_normal_2[i], Bd_FreeSurf_normal_3[i]);
				//cout << "Xi : "<< xi_temp <<"  eta : "<< eta_temp <<" Zeta : "<< zeta_temp <<endl;
				break;
			}

			case TetraAffin:
			{
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((TTetraAffin *)F_K)->SetCell(currentCell);

				((TTetraAffin *)F_K)->GetOuterNormal(JointNumber,xi_1[i],xi_2[i],
											Bd_FreeSurf_normal_1[i], Bd_FreeSurf_normal_2[i], Bd_FreeSurf_normal_3[i]);
				break;
			}
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function get_surface_normals " <<endl;
				exit(0);
			}
		}
		norm1 = Bd_FreeSurf_normal_1[i]; norm2 = Bd_FreeSurf_normal_2[i] ; norm3 = Bd_FreeSurf_normal_3[i];
		normLen = sqrt(norm1*norm1 + norm2*norm2 + norm3*norm3);
		double norm_1  = norm1/normLen; 
		double norm_2 = norm2/normLen; 
		double norm_3 = norm3/normLen; 
		Bd_FreeSurf_normal_1[i] = norm_1; Bd_FreeSurf_normal_2[i] = norm_2 ; Bd_FreeSurf_normal_3[i] = norm_3;

		// ------- NOmal Correction for residual Values
		if ( fabs(Bd_FreeSurf_normal_1[i]) < 1e-7 ) Bd_FreeSurf_normal_1[i] = 0.0;
		if ( fabs(Bd_FreeSurf_normal_2[i]) < 1e-7 ) Bd_FreeSurf_normal_2[i] = 0.0;
		if ( fabs(Bd_FreeSurf_normal_3[i]) < 1e-7 ) Bd_FreeSurf_normal_3[i] = 0.0;


		// cout << " Tangent Vectors ----------- " <<endl;
		// cout << " Cell number : " << freeSurfaceCells[i] <<endl;
		// cout << " Cell Mapping nu : "<< m_h_FreeSurfCellMapping[i] <<endl;


		//////////// -END- CODE BLOCK A2 - Calculate Joint Normal using "MESH" velocity FE Space  ////////////////////
	}

 	// cout << " Size of freeSurfaceCells : " << freeSurfaceCells.size()<<endl;
	// cout << " Size of freeSurfaceJoints: " << freeSurfaceJoints.size()<<endl;
	// cout << " Size of freeSurfaceVertex : " << freeSurfaceVertex.size()<<endl;
}


void FE3D_ALE::get_velocityValues_freeSurface(TFEVectFunct3D* MeshVelo_FEvect,TFEVectFunct3D* Velocity_FEvect )
{
	TBaseCell* currentCell;  TVertex* currentVertex;	
	FE3D FEId, FEId_velocity; TFE3D *ele , *ele_Velocity;
	BF3DRefElements RefElement , RefElement_Velocity;	
	RefTrans3D RefTrans , RefTransVelocity;TRefTrans3D *F_K , *F_K_V;
	TBaseFunct3D *bf;
	TFESpace3D *Velocity_FESpace,*Mesh_FESpace;
	TFEFunction3D*  Velocity_FEFunction[3];
	Mesh_FESpace = MeshVelo_FEvect->GetFESpace3D();
	Velocity_FESpace = Velocity_FEvect->GetFESpace3D();

	for ( int i = 0 ; i < N_freeSurfaceVertex ; i ++)
	{
		int vertex_number = i;
		// cout << " --------------------------- vertex number " <<i <<"---------------------- " <<endl;
		// cout << " Local DOF of the Vertex : " << freeSurfaceVertexLocal[i] ;
		//cout  << " Global DOF of the Vertex : " << freeSurfaceVertex[i] <<endl;
		int cellNr = freeSurfaceCells[i];
		currentCell =  coll->GetCell(cellNr);

		// Parameters for the Mesh Velocity FESPACE	---- /////
		FEId = Mesh_FESpace->GetFE3D(cellNr, currentCell);
		ele = TFEDatabase3D::GetFE3D(FEId);
		RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
		RefTrans3D referenceTransformation =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////


		// Parameters for the  Velocity FESPACE	---- /////
		FEId_velocity = Velocity_FESpace->GetFE3D(cellNr, currentCell);
		ele_Velocity = TFEDatabase3D::GetFE3D(FEId_velocity);
		RefElement_Velocity = TFEDatabase3D::GetRefElementFromFE3D(FEId_velocity);
		RefTrans3D RefTransVelocity =  ele_Velocity->GetRefTransID();
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////

		//----------- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //
		int* GlobalNumbers_Velocity = Velocity_FESpace->GetGlobalNumbers();
		int* BeginIndex_velocity = Velocity_FESpace->GetBeginIndex();
		int* Numbers_Velocity  = GlobalNumbers_Velocity + BeginIndex_velocity[cellNr];

		int* GlobalNumbers_Mesh = Mesh_FESpace->GetGlobalNumbers();
		int* BeginIndex_Mesh = Mesh_FESpace->GetBeginIndex();
		int* Numbers_Mesh = GlobalNumbers_Mesh + BeginIndex_Mesh[cellNr];
		//----END---- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //

		double* u_values = new double[3]();
		double* meshVelocityValues = new double[3]();
		//Get the Global Co-ordinates of the DOF in Velocity FESpace
		TBaseCell* VelocityCell =  MovCells[i];
		int freeSurface_globalID = freeSurfaceVertex[i];
		switch(RefTransVelocity)
		{
			case HexaTrilinear: 
			{
				int*  BeginIndex = Velocity_FEvect->GetFESpace3D()->GetBeginIndex();
				int* GlobalNumbers = Velocity_FEvect->GetFESpace3D()->GetGlobalNumbers();

				FE3D FE_ID = Velocity_FEvect->GetFESpace3D()->GetFE3D(freeSurfaceCells[vertex_number], currentCell);
				TFE3D* FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
				RefTrans = FE_Obj->GetRefTransID();
				TFEDatabase3D::SetCellForRefTrans(currentCell, RefTrans);
				// ------ Calculate the Values of the Velocity Function at the Nodal Points   ------- //
				bf = ele_Velocity->GetBaseFunct3D();
				int N_BaseFunct = bf->GetDimension();

				double* uref = new double[N_BaseFunct];

				// bf->GetDerivatives(D000, xiq, etaq, zetaq, uref);
				bf->GetDerivatives(D000,xi_freeSurfaceNodes[vertex_number],
						eta_freeSurfaceNodes[vertex_number], zeta_freeSurfaceNodes[vertex_number], uref);


				double u = 0,val;
				int* Numbers = GlobalNumbers + BeginIndex[freeSurfaceCells[vertex_number]];
				
				for ( int k = 0 ; k < 3 ; k++)
				{
					double* Values = Velocity_FEvect->GetComponent(k)->GetValues();
					for(int j=0;j<N_BaseFunct;j++)
					{
						val = Values[Numbers[j]];
						// cout << j << " " << val << endl;
						meshVelocityValues[k]  +=  uref[j]*val;
					} 
				}
				delete[] uref;
				break;
			}

			case TetraAffin:
			{
				int*  BeginIndex = Velocity_FEvect->GetFESpace3D()->GetBeginIndex();
				int* GlobalNumbers = Velocity_FEvect->GetFESpace3D()->GetGlobalNumbers();

				FE3D FE_ID = Velocity_FEvect->GetFESpace3D()->GetFE3D(freeSurfaceCells[vertex_number], currentCell);
				TFE3D* FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
				RefTrans = FE_Obj->GetRefTransID();
				TFEDatabase3D::SetCellForRefTrans(currentCell, RefTrans);
				// ------ Calculate the Values of the Velocity Function at the Nodal Points   ------- //
				bf = ele_Velocity->GetBaseFunct3D();
				int N_BaseFunct = bf->GetDimension();

				double* uref = new double[N_BaseFunct];

				// bf->GetDerivatives(D000, xiq, etaq, zetaq, uref);
				bf->GetDerivatives(D000,xi_freeSurfaceNodes[vertex_number],
						eta_freeSurfaceNodes[vertex_number], zeta_freeSurfaceNodes[vertex_number], uref);


				double u = 0,val;
				int* Numbers = GlobalNumbers + BeginIndex[freeSurfaceCells[vertex_number]];
				
				for ( int k = 0 ; k < 3 ; k++)
				{
					double* Values = Velocity_FEvect->GetComponent(k)->GetValues();
					for(int j=0;j<N_BaseFunct;j++)
					{
						val = Values[Numbers[j]];
						// cout << j << " " << val << endl;
						meshVelocityValues[k]  +=  uref[j]*val;
					} 
				}

				delete[] uref;
				break;
			}
			
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function: get_surface_normals " <<endl;
				exit(0);
			}
		}

		// Velocity_FEvect->FindValueLocal(MovCells[i], freeSurfaceCells[i], 
									//   x_coord,y_c		// 
		// double x_coord, y_coord,z_coord;
		// Mesh_FESpace->GetDOFPosition(freeSurfaceVertex[i], x_coord,y_coord, z_coord);
		// double uVal[4] ={0}, vVal[4] = {0} , wVal[4] = {0};
		// Velocity_FEvect->GetComponent(0)->FindGradientLocal(MovCells[i],freeSurfaceCells[i],x_coord, y_coord, z_coord, uVal);
		// Velocity_FEvect->GetComponent(1)->FindGradientLocal(MovCells[i],freeSurfaceCells[i],x_coord, y_coord, z_coord, vVal);
		// Velocity_FEvect->GetComponent(2)->FindGradientLocal(MovCells[i],freeSurfaceCells[i],x_coord, y_coord, z_coord, wVal);
		// int checker = 0; 
		// for ( int k = 0 ; k < 3 ; k++)
		// {
		// 	if ( (fabs( uVal[0] - meshVelocityValues[0]) > 1e-9 )
		// 		|| (fabs( vVal[0] - meshVelocityValues[1]) > 1e-9 )
		// 		|| (fabs( wVal[0] - meshVelocityValues[2]) > 1e-9 ))
		// 	{
		// 		checker++;
		// 	}
		// }
		// if(checker) cout << " ERORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR" <<endl;



		u_values[0] = meshVelocityValues[0];
		u_values[1] = meshVelocityValues[1];
		u_values[2] = meshVelocityValues[2];


		//////////// - END -- CODE BLOCK B - To get the Normal Vectors and Velocity Values at the nodal Points ///////////////////////
		// cout << " Velocity Vector  : [ " << u_values[0] << ", " << u_values[1] << ", " << u_values[2] << " ] " <<endl;
		// cout << "node : " <<  "   Norm Vector  : [ " << Bd_FreeSurf_normal_1[i] << ", " << Bd_FreeSurf_normal_2[i] << ", " << Bd_FreeSurf_normal_3[i] << " ]  "  <<endl;
		double alpha = 0;
		alpha = u_values[0]*Bd_FreeSurf_normal_1[i] + u_values[1]*Bd_FreeSurf_normal_2[i] + u_values[2]*Bd_FreeSurf_normal_3[i];	
		// cout << " Alpha : " << alpha <<endl;

		

		// ASsign the ith component of Mesh Velocity as ( u1.n1 + u2.n2 + u3.n3 ) * Normal
		// HARDCODED - For 
		meshVelocityatFreeSurface_1[vertex_number] = alpha*Bd_FreeSurf_normal_1[i];   // Value of w1 of Global DOF freeSurfaceVertex[i]
		meshVelocityatFreeSurface_2[vertex_number] = alpha*Bd_FreeSurf_normal_2[i];   // Value of w2 of Global DOF freeSurfaceVertex[i]
		meshVelocityatFreeSurface_3[vertex_number] = alpha*Bd_FreeSurf_normal_3[i];   // Value of w3 of Global DOF freeSurfaceVertex[i]
		
		// meshVelocityatFreeSurface_1[vertex_number] = u_values[0];   // Value of w1 of Global DOF freeSurfaceVertex[i]
		// meshVelocityatFreeSurface_2[vertex_number] = u_values[1];   // Value of w2 of Global DOF freeSurfaceVertex[i]
		// meshVelocityatFreeSurface_3[vertex_number] = u_values[2];   // Value of w3 of Global DOF freeSurfaceVertex[i]


		delete[] u_values;
			
		}
}


void FE3D_ALE::get_surface_normals_slipBoundary(TFEVectFunct3D* MeshVelo_FEvect )
{
	//variable Declarations
	TBaseCell* currentCell;  TVertex* currentVertex;	
	FE3D FEId, FEId_velocity; TFE3D *ele;
	BF3DRefElements RefElement;	
	RefTrans3D RefTrans , RefTransVelocity;TRefTrans3D *F_K;
	TBaseFunct3D *bf;
	TFESpace3D *Mesh_FESpace;
	TFEFunction3D*  Mesh_FEFunction[3];


	//  Get the Compunents of MeshVelo_fevect
	int meshVelComp = MeshVelo_FEvect->GetN_Components();

	//get FE Space for Velocity and 
	Mesh_FESpace = MeshVelo_FEvect->GetFESpace3D();	

	double* MeshVelocityArray = MeshVelo_FEvect->GetComponent(0)->GetValues();


	// Get the Lengths of the Value Arrays of Respective FEFunction Values
	int MeshVeloLength = MeshVelo_FEvect->GetComponent(0)->GetLength();

	double norm1 = 0.,norm2 = 0.,norm3 = 0.,normLen = 0.;

	// cout << " N BD CELLS : " << N_bd_FreeSlip_Vertex<<endl;

	//Calculate normal for Each Node in the Joint. 
	// Hard Code the Values of Joint id's for the given HEXAHEADRAL - TRILINEAR
	for ( int i = 0 ; i < N_bd_FreeSlip_Vertex ; i ++)
	{
		int vertex_number = i;
		// cout << " --------------------------- vertex number " <<i <<"---------------------- " <<endl;
		// cout << " Local DOF of the Vertex : " << freeSurfaceVertexLocal[i] ;
		//cout  << " Global DOF of the Vertex : " << freeSurfaceVertex[i] <<endl;
		int cellNr = Bd_FreeSlip_Cells[i];
		currentCell =  coll->GetCell(cellNr);

		// Parameters for the Mesh Velocity FESPACE	---- /////
		FEId = Mesh_FESpace->GetFE3D(cellNr, currentCell);
		ele = TFEDatabase3D::GetFE3D(FEId);
		RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
		RefTrans3D referenceTransformation =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////


		//----------- Get the Local to Global Numbering of the  MESH Fe Space  ----  //
		int* GlobalNumbers_Mesh = Mesh_FESpace->GetGlobalNumbers();
		int* BeginIndex_Mesh = Mesh_FESpace->GetBeginIndex();
		int* Numbers_Mesh = GlobalNumbers_Mesh + BeginIndex_Mesh[cellNr];
		//----END---- Get the Local to Global Numbering of the MESH Fe Space  ----  //


		// SETUP the Values of s and T that needs to be sent to "GetOuterNormal"/ getTangentVectors to get the normal at the point
		double x_coord ,y_coord,z_coord;
		//Mesh_FESpace->GetDOFPosition(freeSurfaceVertex[i], x_coord,y_coord, z_coord);
		// cout << "X : "<< x_coord <<"  Y : "<< y_coord <<"   Z : "<< z_coord <<endl;
		int JointNumber = Bd_FreeSlip_Joints[i];
		
		//////////// CODE BLOCK A2 - Calculate Joint Normal using "MESH" velocity FE Space  ////////////////////
		////// Note : This Code Block is Repeated Block of Code to calculate Normal of a Face at the mesh Face //////
		////// Note : This Code is for "MESH" Velocity FE Space only /////////////////
		double t11 = 0, t12 = 0, t13= 0,t21= 0,t22= 0,t23= 0;
		double xi, eta, zeta;
		double xi_1,xi_2;
		switch(referenceTransformation)     // Reftrans of MESH Velocity
    	{
			case HexaTrilinear: // HEXATRILINEAR 
			{	
				//RefTrans = HexaTrilinear;
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((THexaTrilinear *)F_K)->SetCell(currentCell);
				int localNodeNUmber = Bd_FreeSlip_VertexLocal[vertex_number];
				( (THexaTrilinear *) F_K )->GetRefvaluesfromLocalNodeNumber(FEId,localNodeNUmber,xi,eta,zeta);
				( (THexaTrilinear *) F_K )->GetRefValuesfromJointid(JointNumber,xi,eta,zeta,xi_1,xi_2);
				( (THexaTrilinear *) F_K )->GetTangentVectors(JointNumber,xi_1,xi_2,t11,t12,t13,t21,t22,t23);
				//cout << "Xi : "<< xi_temp <<"  eta : "<< eta_temp <<" Zeta : "<< zeta_temp <<endl;
				break;
			}

			case TetraAffin:
			{
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((TTetraAffin *)F_K)->SetCell(currentCell);
				int localNodeNUmber = Bd_FreeSlip_VertexLocal[vertex_number];
				( (TTetraAffin *) F_K )->GetRefvaluesfromLocalNodeNumber(FEId,localNodeNUmber,xi,eta,zeta);
				
				// The below Step is redundant , as the s and t ( xi_1 and xi_2 ) values will not be used for calculation of 
				// normals at the Tetraheadral Cell, However to maintain code consistency with hexalinear, we will assign the values to any of the 
				// twi reference co-ordinates 
				xi_1 = xi; xi_2 = eta;
				// end

				( (TTetraAffin *) F_K )->GetTangentVectors(JointNumber,xi_1,xi_2,t11,t12,t13,t21,t22,t23);
				break;
			}
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function:  get_surface_normals " <<endl;
				exit(0);
			}
		}

		// Normalise the Tangent Vectors
		double tang1Len = sqrt(t11*t11 + t12*t12 + t13*t13);
		double tang2Len = sqrt(t21*t21 + t22*t22 + t23*t23);
		t11 = t11/tang1Len; t12 = t12/tang1Len ; t13 = t13/tang1Len;
		t21 = t21/tang2Len ; t22 = t22/tang2Len ; t23 = t23/tang2Len;
		

		// save the normal  Values to the tangent Vectors array
		Bd_TangentA_1[i] = t11;Bd_TangentA_2[i] = t12;Bd_TangentA_3[i] = t13;
		Bd_TangentB_1[i] = t21;Bd_TangentB_2[i] = t22;Bd_TangentB_3[i] = t23;

		// NOrmlaise the Normal Vectors
		norm1 = t12*t23 - t13*t22;
  		norm2 = t13*t21 - t11*t23;
  		norm3 = t11*t22 - t12*t21;
		normLen = sqrt(norm1*norm1 + norm2*norm2 + norm3*norm3);
		
		// save the normal  Values to the tangent Vectors array
		Bd_normal_1[i]  = norm1/normLen; 
		Bd_normal_2[i] = norm2/normLen; 
		Bd_normal_3[i] = norm3/normLen; 

		if(fabs(Bd_normal_1[i]) < 1e-7)  Bd_normal_1[i] = 0.0;
		if(fabs(Bd_normal_2[i]) < 1e-7)  Bd_normal_2[i] = 0.0;
		if(fabs(Bd_normal_3[i]) < 1e-7)  Bd_normal_3[i] = 0.0;
		

		// if( (fabs(Bd_normal_1[i]) - 1.0) < 1e-7) 
		// 	{ Bd_normal_1[i] > 0 ?  Bd_normal_1[i] = 1.0  : Bd_normal_1[i] = -1.0  ; }
		// if((fabs(Bd_normal_2[i]) - 1.0) < 1e-7) 
		// 	{ Bd_normal_2[i] > 0 ?  Bd_normal_2[i] = 1.0  : Bd_normal_2[i] = -1.0  ; }
		// if((fabs(Bd_normal_3[i]) - 1.0 )< 1e-7) 
		// 	{ Bd_normal_3[i] > 0 ?  Bd_normal_3[i] = 1.0  : Bd_normal_3[i] = -1.0  ; }

		// cout << " DOF : " << Bd_FreeSlip_Vertex[i] <<endl;
		// cout << " [t1,t2,t3] : "<< t11 <<", " <<t12 <<", "<<t13 <<endl;
		// cout << " [t1,t2,t3] : "<< t21 <<", " <<t22 <<", "<<t23 <<endl;
		// cout <<" DOF : " <<Bd_FreeSlip_Vertex[i] << " [n1,n2,n3] : "<< Bd_normal_1[i] << ", "<<Bd_normal_2[i] <<", "<< Bd_normal_3[i]<<endl;

		// cout <<" Vertex : " << Bd_FreeSlip_Vertex[vertex_number] << " Norm Vector  : [ " << Bd_normal_1[i] << ", " << Bd_normal_2[i] << ", " << Bd_normal_3[i] << " ]  Norm len : "  << normLen <<endl;
		//////////// -END- CODE BLOCK A2 - Calculate Joint Normal using MESH FE Space  ////////////////////
	}
}


void FE3D_ALE::get_surface_normals_slipBoundary_EdgeNodes(TFEVectFunct3D* MeshVelo_FEvect )
{
	//variable Declarations
	TBaseCell* currentCell;  TVertex* currentVertex;	
	FE3D FEId, FEId_velocity; TFE3D *ele;
	BF3DRefElements RefElement;	
	RefTrans3D RefTrans , RefTransVelocity;TRefTrans3D *F_K;
	TBaseFunct3D *bf;
	TFESpace3D *Mesh_FESpace;
	TFEFunction3D*  Mesh_FEFunction[3];

	//  Get the Compunents of MeshVelo_fevect
	int meshVelComp = MeshVelo_FEvect->GetN_Components();

	//get FE Space for Velocity and 
	Mesh_FESpace = MeshVelo_FEvect->GetFESpace3D();	

	double* MeshVelocityArray = MeshVelo_FEvect->GetComponent(0)->GetValues();

	// Get the Lengths of the Value Arrays of Respective FEFunction Values
	int MeshVeloLength = MeshVelo_FEvect->GetComponent(0)->GetLength();

	double norm1 = 0.,norm2 = 0.,norm3 = 0.,normLen = 0.;
	

	// cout << " N BD EDGE CELLS : " << N_bd_EdgeFreeSlip_Vertex<<endl;


	//Calculate normal for Each Node in the Joint. 
	// Hard Code the Values of Joint id's for the given HEXAHEADRAL - TRILINEAR
	for ( int i = 0 ; i < N_bd_EdgeFreeSlip_Vertex ; i ++)
	{
		int vertex_number = i;
		// cout << " --------------------------- vertex number " <<i <<"---------------------- " <<endl;
		// cout << " Local DOF of the Vertex : " << freeSurfaceVertexLocal[i] ;
		//cout  << " Global DOF of the Vertex : " << freeSurfaceVertex[i] <<endl;
		int cellNr = Bd_FreeSlip_Cells[i];
		currentCell =  coll->GetCell(cellNr);

		// Parameters for the Mesh Velocity FESPACE	---- /////
		FEId = Mesh_FESpace->GetFE3D(cellNr, currentCell);
		ele = TFEDatabase3D::GetFE3D(FEId);
		RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);
		RefTrans3D referenceTransformation =  TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);
		// --END-- Parameters for the Mesh Velocity FESPACE	---- /////


		//----------- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //
		int* GlobalNumbers_Mesh = Mesh_FESpace->GetGlobalNumbers();
		int* BeginIndex_Mesh = Mesh_FESpace->GetBeginIndex();
		int* Numbers_Mesh = GlobalNumbers_Mesh + BeginIndex_Mesh[cellNr];
		//----END---- Get the Local to Global Numbering of the Velocity FE Space and MESH Fe Space  ----  //


		// SETUP the Values of s and T that needs to be sent to "GetOuterNormal"/ getTangentVectors to get the normal at the point

		int JointNumber1 = Bd_EdgeFreeSlip_Joints1[i];
		int JointNumber2 = Bd_EdgeFreeSlip_Joints2[i];

		
		//////////// CODE BLOCK A2 - Calculate Joint Normal using "MESH" velocity FE Space  ////////////////////
		////// Note : This Code Block is Repeated Block of Code to calculate Normal of a Face at the mesh Face //////
		////// Note : This Code is for "MESH" Velocity FE Space only /////////////////
		double t11 = 0, t12 = 0, t13= 0,t21= 0,t22= 0,t23= 0;
		double xi, eta, zeta;
		double xi_1,xi_2;  // For Normal 
		double xi_3,xi_4;  // For Normal 
		switch(referenceTransformation)     // Reftrans of MESH Velocity
    	{
			case HexaTrilinear: // HEXATRILINEAR 
			{	
				//RefTrans = HexaTrilinear;
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((THexaTrilinear *)F_K)->SetCell(currentCell);
				int localNodeNUmber = Bd_FreeSlip_VertexLocal[vertex_number];
				( (THexaTrilinear *) F_K )->GetRefvaluesfromLocalNodeNumber(FEId,localNodeNUmber,xi,eta,zeta);
				( (THexaTrilinear *) F_K )->GetNormalVectors(JointNumber1,xi,eta,zeta,Bd_edge_normalA_1[i],Bd_edge_normalA_2[i],
															Bd_edge_normalA_3[i],true);
				( (THexaTrilinear *) F_K )->GetNormalVectors(JointNumber2,xi,eta,zeta,Bd_edge_normalB_1[i],Bd_edge_normalB_2[i],
															Bd_edge_normalB_3[i],true);
				( (THexaTrilinear *) F_K )->GetTangentVectors(JointNumber1,xi_1,xi_2,t11,t12,t13,t21,t22,t23);
				//cout << "Xi : "<< xi_temp <<"  eta : "<< eta_temp <<" Zeta : "<< zeta_temp <<endl;
				break;
			}

			case TetraAffin:
			{
				F_K = TFEDatabase3D::GetRefTrans3D(referenceTransformation);
				((TTetraAffin *)F_K)->SetCell(currentCell);
				int localNodeNUmber = Bd_FreeSlip_VertexLocal[vertex_number];
				((TTetraAffin *)F_K)->GetOuterNormal(JointNumber1,xi,eta,Bd_edge_normalA_1[i],Bd_edge_normalA_2[i],
															Bd_edge_normalA_3[i]);

				((TTetraAffin *)F_K)->GetOuterNormal(JointNumber2,xi,eta,Bd_edge_normalB_1[i],Bd_edge_normalB_2[i],
															Bd_edge_normalB_3[i]);

				( (THexaTrilinear *) F_K )->GetTangentVectors(JointNumber1,xi_1,xi_2,t11,t12,t13,t21,t22,t23);
				break;	
			}
			default:
			{
				cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
				cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
				cout << " ERROR LOCATION : Class DeformMesh3D , Function get_surface_normals " <<endl;
				exit(0);
			}
		}

		// Normalise the Tangent Vectors
		double tang1Len = sqrt(t11*t11 + t12*t12 + t13*t13);
		t11 = t11/tang1Len; t12 = t12/tang1Len ; t13 = t13/tang1Len;

		
		// save the normal  Values to the tangent Vectors array
		Bd_edge_TangentA_1[i] = t11;Bd_edge_TangentA_2[i] = t12;Bd_edge_TangentA_3[i] = t13;


		// if( fabs(Bd_edge_normalA_1[i]) - 1.0 < 1e-7) { Bd_edge_normalA_1[i] > 0 ?  Bd_edge_normalA_1[i] = 1.0  : Bd_edge_normalA_1[i] = -1.0  ; }
		// if(fabs(Bd_edge_normalA_2[i]) - 1.0 < 1e-7) { Bd_edge_normalA_2[i] > 0 ?  Bd_edge_normalA_2[i] = 1.0  : Bd_edge_normalA_2[i] = -1.0  ; }
		// if(fabs(Bd_edge_normalA_3[i]) - 1.0 < 1e-7) { Bd_edge_normalA_3[i] > 0 ?  Bd_edge_normalA_3[i] = 1.0  : Bd_edge_normalA_3[i] = -1.0  ; }

		// if(fabs(Bd_edge_normalB_1[i]) - 1.0 < 1e-7) { Bd_edge_normalB_1[i] > 0 ?  Bd_edge_normalB_1[i] = 1.0  : Bd_edge_normalB_1[i] = -1.0  ; }
		// if(fabs(Bd_edge_normalB_2[i]) - 1.0 < 1e-7) { Bd_edge_normalB_2[i] > 0 ?  Bd_edge_normalB_2[i] = 1.0  : Bd_edge_normalB_2[i] = -1.0  ; }
		// if(fabs(Bd_edge_normalB_3[i]) - 1.0 < 1e-7) { Bd_edge_normalB_3[i] > 0 ?  Bd_edge_normalB_3[i] = 1.0  : Bd_edge_normalB_3[i] = -1.0  ; }

		// cout << " Norm Vector  : [ " << Bd_normal_1[i] << ", " << Bd_normal_2[i] << ", " << Bd_normal_3[i] << " ]  Norm len : "  << normLen <<endl;
		//////////// -END- CODE BLOCK A2 - Calculate Joint Normal using MESH FE Space  ////////////////////
	}
}


//  ---------    ALE  GRID TO DATA ---------------------------- //
void FE3D_ALE::move_mesh_ale(TFEVectFunct3D* MeshVelocityVectFunction3D,double* meshVelocity_old , double time_step)
{
	// Create a Vect Fucntion 2D to get the Mesh Co ordinates 
	TFESpace3D* Gridfespace = MeshVelocityVectFunction3D->GetFESpace3D();
	int N_DOF = Gridfespace->GetN_DegreesOfFreedom();
	double* gridpos = new double[3*N_DOF]();
	TFEVectFunct3D* gridFEVectfunc3D = new TFEVectFunct3D(Gridfespace, (char*)"C", (char*)"C", gridpos, N_DOF, 3);

	// value array of current Mesh velocity
	double* meshVelocityvalues = MeshVelocityVectFunction3D->GetValues();

	// get the Co-ordinates as data array using grid to data
	gridFEVectfunc3D->GridToData();

	double maxDispValue = 0;
	double dispValue;
	double avgDispValue = 0; 
	
	// Now the increment in New position of each Vertes is given by ( (W_new + W_old * 0.5) * time_step )
	for ( int i = 0 ; i < 3*N_DOF ; i ++ )
	{
		dispValue   	=  (meshVelocityvalues[i] + meshVelocity_old[i])*0.5 * time_step;
		avgDispValue   +=  (meshVelocityvalues[i] + meshVelocity_old[i])*0.5 * time_step;
		if (dispValue > maxDispValue)  maxDispValue = dispValue;
		gridpos[i] += (meshVelocityvalues[i] + meshVelocity_old[i])*0.5 * time_step;
	}

	avgDispValue   /= (3 * N_DOF) ; 
	// for ( int k = 0 ; k < 3*N_DOF ; k++)
	// 	cout << "  Grid pos ["<<k<<"] : "  << gridpos[k]  <<endl;
	
	// ------------ MESH QUALITY PARAMETER _ ESTIMATIONS ----------------------------- //

	double minShortEdge = coll->GetShortestEdgeinCollection();
	double minVolume    = coll->GetminVolumeCollection();

	double minShortedgeRelative  =       ( 1 -  (InitialMeshMinDiameter - minShortEdge)/ InitialMeshMinDiameter)  * 100 ;

	double minVolumeRelative      =     ( 1 -  (InitialMeshMinVolume - minVolume)/ InitialMeshMinVolume)  * 100 ;

	if( fabs(minShortEdge) < 1e-7 ){
	cout << " ERROR : Mesh has collapsed while moving " <<endl;
	cout << " INFO  : Class - FE3D_ALE ; Function - movemeshale " <<endl;
	exit(0);   // THIVIN - EXIT
	}
	
	if(fabs(minVolume) < 1e-7 ){
	cout << " ERROR : Mesh has collapsed while moving " <<endl;
	cout << " INFO  : Class - FE3D_ALE ; Function - movemeshale " <<endl;
	exit(0);   // THIVIN - EXIT
	}

	// ------------ MESH QUALITY PARAMETER _ ESTIMATIONS ----------------------------- //


	// Now Move the Value from the array to Grid ( vertex Co ordinates )
	gridFEVectfunc3D->DataToGrid();
	
	
	delete[] gridpos;
	
	
	
	cout << " MESH MOVEMENT -- Max_disp :  " << setw(14)<< maxDispValue <<setw(14)<< avgDispValue <<endl;
	cout << " MESH MOVEMENT -- MInShort Edge : " << setw(14) << minShortEdge <<  setw(14) << minShortedgeRelative <<endl;
	cout << " MESH MOVEMENT -- MIn Vol Rel : "<< setw(14) << minVolume << setw(14)<< minVolumeRelative <<endl;


	// cout << " CLASS : DEFORM MESH 3D --  ' Mesh Has been Moved ' " <<endl;

}

// Impose the Free Slip Boundary Condition during Mesh Movement
// THis Function should be called before the Direct Solver Part
void FE3D_ALE::impose_FreeSlip_BoundaryCondition(TSquareStructure3D *sqstructure,double* a11, double* a12,double* a13,
												double* a21, double* a22, double* a23,
												double* a31, double* a32, double* a33 , double* rhs,int length, int N_ActiveBound)
{
	int* RowPtr = sqstructure->GetRowPtr();
	int* ColPtr = sqstructure->GetKCol();

	double Row1_val1 = 0,Row1_val2 = 0,Row1_val3 = 0,Row2_val1 = 0,
		   Row2_val2 = 0, Row2_val3 = 0,Row3_val1 = 0,Row3_val2 = 0,Row3_val3 = 0;
	

	
	for ( int vertexNo = 0 ; vertexNo < N_bd_FreeSlip_Vertex ; vertexNo++)
	{
		// Replace the 1st row with the Normal Values Directly
		int DOF   = Bd_FreeSlip_Vertex[vertexNo];
		int Begin = RowPtr[DOF];
		int End   =  RowPtr[DOF + 1];

		// cout << " DOF : "<< DOF <<endl;
		
		for ( int rowIndex = Begin ; rowIndex < End ; rowIndex ++ )
		{
			// First Row of Blocks
			// Check for the Condition , if any normal is zero or close to zero, then we need to update the Values in the respective 
			// blocks

			// Update in 1st Block
			if( fabs( Bd_normal_1[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a11[rowIndex] = Bd_normal_1[vertexNo];
					a12[rowIndex] = Bd_normal_2[vertexNo];
					a13[rowIndex] = Bd_normal_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a11[rowIndex] = 0;
					a12[rowIndex] = 0;
					a13[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}	

				rhs[DOF + 0*length] = 0.0;
			}

			// Update in 2nd Block
			else if( fabs( Bd_normal_2[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a21[rowIndex] = Bd_normal_1[vertexNo];
					a22[rowIndex] = Bd_normal_2[vertexNo];
					a23[rowIndex] = Bd_normal_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a21[rowIndex] = 0;
					a22[rowIndex] = 0;
					a23[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}	
				
				rhs[DOF + 1*length] = 0.0;
			}

			// Update in 3rd Block
			else if( fabs( Bd_normal_3[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a31[rowIndex] = Bd_normal_1[vertexNo];
					a32[rowIndex] = Bd_normal_2[vertexNo];
					a33[rowIndex] = Bd_normal_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a31[rowIndex] = 0;
					a32[rowIndex] = 0;
					a33[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}			
				rhs[DOF + 2*length] = 0.0;
			}
			
			else
			{
				cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl;
				cout << " Error in Normal Values of FreeSlip DOF "<< DOF << " ( No normal component value is greater than 1 ) -vert NUm :" <<vertexNo <<endl;
				cout << " Error in File FE3D_ALE.C - Function : impose_FreeSlip_BoundaryCondition"<<endl;
				exit(0);
				// THIVIN - EXIT statement
			}
		}

	}

	// ---------------------------------- FOR EDGE DOF ------------------------------------------------ // 
	// It is assumed that the two normals picked for the edge DOF does not have a same non zero component 
	for ( int vertexNo = 0 ; vertexNo < N_bd_EdgeFreeSlip_Vertex ; vertexNo++)
	{
		// Replace the 1st row with the Normal Values Directly
		int DOF 	= 	Bd_EdgeFreeSlip_Vertex[vertexNo];
		int Begin 	= 	RowPtr[DOF];
		int End   	=  	RowPtr[DOF + 1];
		
		for ( int rowIndex = Begin ; rowIndex < End ; rowIndex++ )
		{

			if( fabs( Bd_edge_normalA_1[vertexNo] ) > (1e-2) )
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a11[rowIndex] = Bd_edge_normalA_1[vertexNo];
					a12[rowIndex] = Bd_edge_normalA_2[vertexNo];
					a13[rowIndex] = Bd_edge_normalA_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a11[rowIndex] = 0;
					a12[rowIndex] = 0;
					a13[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}	

				rhs[DOF + 0*length] = 0.0;
			}

			// Update in 2nd Block
			else if( fabs( Bd_edge_normalA_2[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a21[rowIndex] = Bd_edge_normalA_1[vertexNo];
					a22[rowIndex] = Bd_edge_normalA_2[vertexNo];
					a23[rowIndex] = Bd_edge_normalA_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a21[rowIndex] = 0;
					a22[rowIndex] = 0;
					a23[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}	
				
				rhs[DOF + 1*length] = 0.0;
			}

			// Update in 3rd Block
			else if( fabs( Bd_edge_normalA_3[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a31[rowIndex] = Bd_edge_normalA_1[vertexNo];
					a32[rowIndex] = Bd_edge_normalA_2[vertexNo];
					a33[rowIndex] = Bd_edge_normalA_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a31[rowIndex] = 0;
					a32[rowIndex] = 0;
					a33[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}			
				rhs[DOF + 2*length] = 0.0;
			}
			
			else
			{
				cout << " Error in Normal Values of FreeSlip DOF "<< DOF << " ( No normal component value is greater than 1 )" <<endl;
				cout << " Error in File FE3D_ALE.C - Function : impose_FreeSlip_BoundaryCondition"<<endl;
				exit(0);
				// THIVIN - EXIT statement
			}

			if( fabs( Bd_edge_normalB_2[vertexNo] ) > (1e-2) )
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a21[rowIndex] = Bd_edge_normalB_1[vertexNo];
					a22[rowIndex] = Bd_edge_normalB_2[vertexNo];
					a23[rowIndex] = Bd_edge_normalB_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a21[rowIndex] = 0;
					a22[rowIndex] = 0;
					a23[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}	
				rhs[DOF + 1*length] = 0.0;
			}

			// Update in 2nd Block
			else if( fabs( Bd_edge_normalB_1[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a11[rowIndex] = Bd_edge_normalB_1[vertexNo];
					a12[rowIndex] = Bd_edge_normalB_2[vertexNo];
					a13[rowIndex] = Bd_edge_normalB_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a11[rowIndex] = 0;
					a12[rowIndex] = 0;
					a13[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}	
				
				rhs[DOF + 0*length] = 0.0;
			}

			// Update in 3rd Block
			else if( fabs( Bd_edge_normalB_3[vertexNo] ) > (1e-2))
			{
				if(ColPtr[rowIndex] == DOF)
				{
					a31[rowIndex] = Bd_edge_normalB_1[vertexNo];
					a32[rowIndex] = Bd_edge_normalB_2[vertexNo];
					a33[rowIndex] = Bd_edge_normalB_3[vertexNo];
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}
				else
				{
					a31[rowIndex] = 0;
					a32[rowIndex] = 0;
					a33[rowIndex] = 0;
					// cout << a31[rowIndex] << "\t" << a32[rowIndex] << "\t" << a33[rowIndex] << endl;
				}			
				rhs[DOF + 2*length] = 0.0;
			}
			
			else
			{
				cout << " Error in Normal Values of FreeSlip DOF "<< DOF << " ( No normal component value is greater than 1 )" <<endl;
				cout << " Error in File FE3D_ALE.C - Function : impose_FreeSlip_BoundaryCondition"<<endl;
				exit(0);
				// THIVIN - EXIT statement
			}

		}

	}
	
}

void FE3D_ALE::constructGlobalStiffness(TSquareMatrix3D **sqmatrices, int*& RowPtr,int*& KCol, double*& Entries ,
										int&  N_Row,int& N_NNZ, int n_row,int n_column)
{
		int N_Active;
    int N_Rows;
    int N_Entries, begin, end, pos, l;
    int  *rowptr;
    int   *kcol;
    double  *entries, value;
    double t1, t2;

    N_Active = sqmatrices[0]->GetActiveBound();
    N_Rows = sqmatrices[0]->GetN_Rows();

    N_Entries = n_row * n_column * sqmatrices[0]->GetN_Entries();

	#ifdef _CUDA
	
	if(TDatabase::ParamDB->CUDASOLVERFLAG == 1)
	{
		cudaHostAlloc((void**)&Entries, sizeof(double)*N_Entries, cudaHostAllocDefault) ;
    	cudaHostAlloc((void**)&RowPtr, sizeof(int)*((N_Rows * n_row) + 1), cudaHostAllocDefault);
    	cudaHostAlloc((void**)&KCol, sizeof(int)*N_Entries, cudaHostAllocDefault);

	}
	else
	{
	#endif //_CUDA
		
		Entries = new double[N_Entries];
    	RowPtr = new int[N_Rows * n_row + 1];
    	KCol = new int[N_Entries];
	
	#ifdef _CUDA
	}
	#endif //_CUDA
	


    N_Row = N_Rows * n_row;

    pos = 0;
    RowPtr[0] = 0;

    for (int i = 0; i < n_row; ++i)
    {
      for (int row = 0; row < N_Rows; ++row)
      {
        for (int j = 0; j < n_column; ++j)
        {
          // if ( i != j && row >= N_Active ) continue;

          rowptr = sqmatrices[i * n_column + j]->GetRowPtr();
          kcol = sqmatrices[i * n_column + j]->GetKCol();
          entries = sqmatrices[i * n_column + j]->GetEntries();

          begin = rowptr[row];
          end = rowptr[row + 1];

          for (int loc_pos = begin; loc_pos < end; ++loc_pos)
          {
            Entries[pos] = entries[loc_pos];
            KCol[pos] = kcol[loc_pos] + j * N_Rows;
            ++pos;
          }
        }
        RowPtr[i * N_Rows + row + 1] = pos;
      }
    }

	N_NNZ = pos;
}

void FE3D_ALE::AssembleMeshMatrix(TFESpace3D* fespace, TFEVectFunct3D* MeshVelocityVectFunction3D, int FactoriseFlag)
{
	PUSH_RANGE("Mesh Matrix Assemble", 5);
	char UString[] = "T";
	char NameString[] = "name";
	char CString[] = "C";
	TFESpace3D* Meshfespace = MeshVelocityVectFunction3D->GetFESpace3D();
	int N_DOF = Meshfespace->GetN_DegreesOfFreedom();
	int N_Active = Meshfespace->GetActiveBound();
    
    // cout  << "N_Cells = " << coll -> GetN_Cells() << endl;
	// cout << "Degrees of Freedom = " << N_DOF  << "    N_Active = " << N_Active << endl;
	
	double* meshVelocityValues;
	meshVelocityValues = MeshVelocityVectFunction3D->GetValues();

	GlobRhs = new double[3*N_DOF]();
    
    TAuxParam3D *Meshaux = new TAuxParam3D(1, 0, 0, 0, &Meshfespace, NULL, NULL, NULL, NULL, 0, NULL);
    
    int N_Terms = 4;    
  	int *SpacesNumbers = new int[N_Terms](); 
  	int N_Matrices = 9;    	
	int N_RHS = 3;   
	int *rowspace = new int[N_Matrices]();
	int *columnspace = new int[N_Matrices](); 
	int *rhsspace = new int[N_RHS]();
    
    
    MultiIndex3D AllDerivatives[4] = {D000, D100, D010,D001};
     
    TDiscreteForm3D* discreteform = new TDiscreteForm3D(UString, UString, N_Terms, AllDerivatives,
                                        SpacesNumbers, N_Matrices, N_RHS, rowspace, columnspace, rhsspace,
										                    Assembly_poisson_3D, GridCoeffs, NULL); 
    
    // --------------------- START OF MATRIX STRUCTURE DECLARATION -------------------//



	SQMATRICES_GRID[0] = SqmatrixG11;  SQMATRICES_GRID[1] = SqmatrixG12; SQMATRICES_GRID[2] = SqmatrixG13;
	SQMATRICES_GRID[3] = SqmatrixG21;  SQMATRICES_GRID[4] = SqmatrixG22; SQMATRICES_GRID[5] = SqmatrixG23;
	SQMATRICES_GRID[6] = SqmatrixG31;  SQMATRICES_GRID[7] = SqmatrixG32; SQMATRICES_GRID[8] = SqmatrixG33;

	RHS[0] = GlobRhs; RHS[1] = GlobRhs + N_DOF; RHS[2] = GlobRhs + 2*N_DOF;

	Entries[0] = SqmatrixG11->GetEntries(); Entries[1] = SqmatrixG12->GetEntries(); Entries[2] = SqmatrixG13->GetEntries();
	Entries[3] = SqmatrixG21->GetEntries(); Entries[4] = SqmatrixG22->GetEntries(); Entries[5] = SqmatrixG23->GetEntries();
	Entries[6] = SqmatrixG31->GetEntries(); Entries[7] = SqmatrixG32->GetEntries(); Entries[8] = SqmatrixG33->GetEntries();
	
	GridKCol = sqstructure->GetKCol(); GridRowPtr = sqstructure->GetRowPtr();
	
	// ---------------- START OF ASSEMBLY 3D FUNCTION -----------//
	fesp[0] = Meshfespace;   // Type of FE Space to be used for Blocks in A Matrix
    ferhs[0] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
    ferhs[1] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	ferhs[2] = Meshfespace;  // Type of FE Space to be used for Blocks in Rhs Matrix
	
    TAssembleMat3D *MeshMatAssemble = new TAssembleMat3D(1, &Meshfespace, 9, SQMATRICES_GRID, 0, NULL,
														 3, RHS, ferhs, discreteform, GridBoundaryConditions,
														  GridBoundValues, Meshaux);
	MeshMatAssemble->Init();
	MeshMatAssemble->Reset();
	MeshMatAssemble->Assemble3D();

	// N_Active = Non Dirichlet DOF's
	int N_BDDof = N_DOF - N_Active;

	int N_ActiveBoundary = SqmatrixG11->GetActiveBound();
	

	const int N_DOF_perBlock =  SqmatrixG11->GetN_Rows();
	const int N_Active_perBlock = SqmatrixG11->GetActiveBound();
	const int N_Dirichlet_perBlock = N_DOF - N_Active;
	GridRowPtr = sqstructure->GetRowPtr();

	
	// // Memset the Antidiagonal arrays to be zero. 
	memset(SqmatrixG12->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] ) *sizeof(double));
	memset(SqmatrixG13->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG21->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG23->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG31->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));
	memset(SqmatrixG32->GetEntries() + GridRowPtr[N_Active_perBlock],0,( GridRowPtr[N_DOF]-GridRowPtr[N_Active] )*sizeof(double));

	
	// Synchronise the streams which has calculated the normals
	#ifdef _CUDA
	checkCudaErrors(cudaStreamSynchronize(FE3d_stream));
	#endif  // _CUDA
	// Assign the u.n boundary Condition for the Free Slip Slides
	impose_FreeSlip_BoundaryCondition(sqstructure,SqmatrixG11->GetEntries(),SqmatrixG12->GetEntries(),SqmatrixG13->GetEntries(),SqmatrixG21->GetEntries(),SqmatrixG22->GetEntries()
									,SqmatrixG23->GetEntries(),SqmatrixG31->GetEntries(),SqmatrixG32->GetEntries(),SqmatrixG33->GetEntries(),GlobRhs,N_DOF,N_ActiveBoundary);
	POP_RANGE;


	#ifdef _CUDA
	// If Cuda Solver is Selected , Solve here, Else Solve in THe Solve MEsh Matrix Routine
	if(TDatabase::ParamDB->CUDASOLVERFLAG && FactoriseFlag)
	{	
		FactoriseGlobalMeshMatrix(SQMATRICES_GRID, 3, 3);
	}
	#endif  // _CUDA

}

#ifdef _CUDA

void FE3D_ALE::FactoriseGlobalMeshMatrix(TSquareMatrix3D **sqmatrices, int n_row , int n_column)
{
	// Do Refactoring Steps , if the Iteration is not the first Iteration and the Refactorisation count is less than the
	// Value defind in Dat file
	if(N_refactorStepsDone && N_refactorStepsDone <= TDatabase::ParamDB->N_CUDA_REFACTOR_STEPS )
	{

		int* RowPtr;	
		int* KCol;
		double* Entries;
		int N_Row;
		int N_NNZ;
		constructGlobalStiffness(sqmatrices,RowPtr,KCol,Entries,N_Row,N_NNZ,3,3);
		auto start = omp_get_wtime();
		if(0 == (strcmp(TDatabase::ParamDB->CUDASOLVERTYPE, "REFACTORLU")))
		{
    		cudaSolverReFactLU->cudaRefactorize(RowPtr, KCol,Entries,N_Row,N_refactorStepsDone );
		}
		auto end = omp_get_wtime();
		cout << " Refactor time Taken : " << end - start << " seconds " <<endl;

		N_refactorStepsDone++;
		if(TDatabase::ParamDB->CUDASOLVERFLAG == 1)
		{
			checkCudaErrors(cudaFreeHost(Entries));
			checkCudaErrors(cudaFreeHost(KCol));
			checkCudaErrors(cudaFreeHost(RowPtr));
		}
		else
		{
			if(RowPtr) delete[] RowPtr;
			if(KCol) delete[] KCol;
			if(Entries) delete[] Entries;
		}
	}
}
#endif  // _CUDA


#ifdef _CUDA
void FE3D_ALE::SolveMeshMatrix(TFESpace3D* fespace, TFEVectFunct3D* MeshVelocityVectFunction3D)
{
	TFESpace3D* Meshfespace = MeshVelocityVectFunction3D->GetFESpace3D();
	int N_DOF = Meshfespace->GetN_DegreesOfFreedom();
	double *sol = new double[3*N_DOF]();
	// double *rhs =  new double[3*N_DOF]();

	// Transfer the Solution Calculated from "U" as a Dirichlet boundary to the Mesh Velocity Solution
	int cnt = 0;
	int cnt2 = 0;
	for ( int i = 0 ; i < N_freeSurfaceVertex ; i++)
	{	
		int globDOF 	= 	freeSurfaceVertex[i];
		std::vector<int>::iterator edgeDOF_it 	=   std::find(FreeSurfaceEdgeDOFs.begin(), FreeSurfaceEdgeDOFs.end(),globDOF ); 
		if ( edgeDOF_it != FreeSurfaceEdgeDOFs.end())  // Edge DOF
		{
			
			if(fabs(Bd_FreeSurf_normal_2[i]) > 1e-3)
			{
				cnt++;

				// Apply Contraint on X 
				meshVelocityatFreeSurface_2[i] 	+=  meshVelocityatFreeSurface_3[i]*(Bd_FreeSurf_normal_3[i]/Bd_FreeSurf_normal_2[i]) ;
			

				GlobRhs[globDOF           ] =  meshVelocityatFreeSurface_1[i];
				GlobRhs[globDOF +  1*N_DOF] =  meshVelocityatFreeSurface_2[i];
				GlobRhs[globDOF +  2*N_DOF] =  0;
			}
		}

		else
		{
			if(fabs(Bd_FreeSurf_normal_2[i]) > 1e-3)
			{
				cnt2++;
				meshVelocityatFreeSurface_2[i] 	+=  meshVelocityatFreeSurface_3[i]*(Bd_FreeSurf_normal_3[i]/Bd_FreeSurf_normal_2[i]) ;
			

				GlobRhs[globDOF           ] =  meshVelocityatFreeSurface_1[i];
				GlobRhs[globDOF +  1*N_DOF] =  meshVelocityatFreeSurface_2[i];
				GlobRhs[globDOF +  2*N_DOF] =  0;
			}
		}
		

	}

	// cout<<" COUNTERRRRRRR : " << cnt <<endl;
	// cout<<" COUNTERRRRRRR : " << cnt2 <<endl;


	// USe Direct Solver
	if(TDatabase::ParamDB->CUDASOLVERFLAG == 0)
	{
		auto start = high_resolution_clock::now(); 

		PardisoDirectSolver_without_removing_dirichlet(SQMATRICES_GRID, 3, 3, sol, GlobRhs);
		
		auto end = high_resolution_clock::now(); 

		auto duration = duration_cast<milliseconds>(end - start); 
	
		// cout << " Mesh Velocity Solution Norm : " <<Ddot(3*N_DOF,sol,sol)<<endl;
		// cout << " Mesh Velocity Solution Time(Direct) : " <<duration.count() << " ms" <<endl;


		memcpy( MeshVelocityVectFunction3D->GetValues(), sol , SizeOfDouble*3*N_DOF);

		// Release the Matrix Storage Parameters
		for ( int i = 0 ; i < 9 ; i++)
			SQMATRICES_GRID[i]->Reset();

		for (int i_rhs = 0; i_rhs < 3*N_DOF; i_rhs++)
			GlobRhs[i_rhs] = 0;

	}

	else
	{
		// NO refactor Steps are Done , Need to Do LU Factorisattion in CPU
		if(!N_refactorStepsDone)    
		{
			if(0 == (strcmp(TDatabase::ParamDB->CUDASOLVERTYPE, "REFACTORLU")))
			{
				cudaSolverReFactLU->initialiseCudaHandles();

				int* RowPtr;
				int* KCol;
				double* Entries;
				int N_Row;
				int N_NNZ;
				char* reorder = "metis";
				
				constructGlobalStiffness(SQMATRICES_GRID,RowPtr,KCol,Entries,N_Row,N_NNZ,3,3);

				cudaSolverReFactLU->LU_DecompositionHost(RowPtr,KCol,Entries,N_Row,N_NNZ,GlobRhs,reorder);

				cudaStreamSynchronize(cudaSolverReFactLU->stream);

				memcpy(MeshVelocityVectFunction3D->GetValues(),GlobRhs, SizeOfDouble*3*N_DOF);

				cout << " Mesh Velocity Norm --  HOST Factor : " <<Ddot(3*N_DOF,GlobRhs,GlobRhs)<<endl;
				
				N_refactorStepsDone = 1;

				if(TDatabase::ParamDB->CUDASOLVERFLAG == 1)
				{
					
				}
				else
				{
					if(RowPtr)  delete[] RowPtr;
					if(KCol)    delete[] KCol;
					if(Entries) delete[] Entries;
				}
				
			}

		}

		else 
		{
			auto start = omp_get_wtime();
			cudaSolverReFactLU->cudaRefactorSolve(GlobRhs);
			auto end = omp_get_wtime();
			cout << " Solve time : " << end - start << " Seconds " <<endl;

			
			memcpy(MeshVelocityVectFunction3D->GetValues(),GlobRhs, SizeOfDouble*3*N_DOF);

			// cout << " Mesh Velocity Solution Norm LU REFACTOR Factor : " <<Ddot(3*N_DOF,GlobRhs,GlobRhs)<<endl;


			if(N_refactorStepsDone == TDatabase::ParamDB->N_CUDA_REFACTOR_STEPS)
			{
				cudaSolverReFactLU->resetCudaRF();
				N_refactorStepsDone = 0;
				cout << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<endl;
				cout << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<endl;

			}


		}
	}
	
	


	if(GlobRhs) delete[] GlobRhs;
	if(sol) delete[] sol;
}

#endif //_CUDA

double FE3D_ALE::getMaximumElevation(TFEVectFunct3D* MeshVelocityVectFunction3D)
{
	int size = freeSurfaceVertex.size();
	double Ymax = -1 * 1e100;
	for ( int vertNo = 0 ; vertNo < size; vertNo++)
	{
		TBaseCell* currentCell = coll->GetCell(freeSurfaceCells[vertNo]);
		FE3D FE_ID = MeshVelocityVectFunction3D->GetFESpace3D()->GetFE3D(freeSurfaceCells[vertNo], currentCell);
		TFE3D* FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
		RefTrans3D RefTrans = FE_Obj->GetRefTransID();
		TFEDatabase3D::SetCellForRefTrans(currentCell, RefTrans);

		TRefTrans3D* rt = TFEDatabase3D::GetRefTrans3D(RefTrans);

		double x, y , z ;
		double absdetjk[1];

		double *xi, *eta, *zeta;
		xi = &xi_freeSurfaceNodes[vertNo];
		eta = &eta_freeSurfaceNodes[vertNo];
		zeta = &zeta_freeSurfaceNodes[vertNo];

		


		switch(RefTrans)
		{
			case TetraAffin:
			((TTetraAffin *)rt)->SetCell(currentCell);
			((TTetraAffin *)rt)->GetOrigFromRef(1,xi, eta, zeta,
													&x, &y, &z, absdetjk);
			break;
	
			case TetraIsoparametric:
			((TTetraIsoparametric *)rt)->SetCell(currentCell);
			((TTetraIsoparametric *)rt)->GetOrigFromRef(1,xi, eta, zeta,
													&x, &y, &z, absdetjk);
			break;
	
			case HexaAffin:
			((THexaAffin *)rt)->SetCell(currentCell);
			((THexaAffin *)rt)->GetOrigFromRef(1,xi, eta, zeta,
													&x, &y, &z, absdetjk);
			break;
	
			case HexaTrilinear:
			((THexaTrilinear *)rt)->SetCell(currentCell);
			((THexaTrilinear *)rt)->GetOrigFromRef(1,xi, eta, zeta,
													&x, &y, &z, absdetjk);
			break;
	
			default:
			{
				cout << " ERROR " <<endl;
			}
		}
		
		if(fabs(x - 10) < 1e-7 && fabs(z - 1.61) < 0.15)
		{
			Ymax = y;
		}
	
		if (y > Ymax) Ymax = y;

	}

	cout << " TIME STEP :  " << TDatabase::TimeDB->CURRENTTIME << " Max Displacement : " << Ymax <<endl;
	return Ymax;

}


#ifdef _CUDA


void FE3D_ALE::updatevertexCoordinates_freeSurf()
{
  for ( int i = 0 ; i < uniqueCellNumbers.size(); i++)
  {
    TBaseCell* Cell = coll->GetCell(uniqueCellNumbers[i]);

    double* x = new double[8]();
    double* y = new double[8]();
    double* z = new double[8]();
    Cell->GetVertex(0)->GetCoords(x[0], y[0], z[0]);
    Cell->GetVertex(1)->GetCoords(x[1], y[1], z[1]);
    Cell->GetVertex(2)->GetCoords(x[2], y[2], z[2]);
    Cell->GetVertex(3)->GetCoords(x[3], y[3], z[3]); 
    Cell->GetVertex(4)->GetCoords(x[4], y[4], z[4]);
    Cell->GetVertex(5)->GetCoords(x[5], y[5], z[5]);
    Cell->GetVertex(6)->GetCoords(x[6], y[6], z[6]);
    Cell->GetVertex(7)->GetCoords(x[7], y[7], z[7]);

    int displacement = i * m_h_N_VertexPerCell;
    
    checkCudaErrors(cudaMemcpyAsync(X_cord + displacement   ,x, sizeof(double)*m_h_N_VertexPerCell,cudaMemcpyHostToDevice,FE3d_stream));
    checkCudaErrors(cudaMemcpyAsync(Y_cord + displacement   ,y, sizeof(double)*m_h_N_VertexPerCell,cudaMemcpyHostToDevice,FE3d_stream));
    checkCudaErrors(cudaMemcpyAsync(Z_cord + displacement   ,z, sizeof(double)*m_h_N_VertexPerCell,cudaMemcpyHostToDevice,FE3d_stream));
    // cudaDeviceSynchronize();
    delete[] x; delete[] y; delete[] z;
  }
}


void FE3D_ALE::updatevertexCoordinates_freeSlip()
{
	
  for ( int i = 0 ; i < uniqueCellNumbers_freeSlip.size(); i++)
  {
    TBaseCell* Cell = coll->GetCell(uniqueCellNumbers_freeSlip[i]);

    double* x = new double[8]();
    double* y = new double[8]();
    double* z = new double[8]();
    Cell->GetVertex(0)->GetCoords(x[0], y[0], z[0]);
    Cell->GetVertex(1)->GetCoords(x[1], y[1], z[1]);
    Cell->GetVertex(2)->GetCoords(x[2], y[2], z[2]);
    Cell->GetVertex(3)->GetCoords(x[3], y[3], z[3]); 
    Cell->GetVertex(4)->GetCoords(x[4], y[4], z[4]);
    Cell->GetVertex(5)->GetCoords(x[5], y[5], z[5]);
    Cell->GetVertex(6)->GetCoords(x[6], y[6], z[6]);
    Cell->GetVertex(7)->GetCoords(x[7], y[7], z[7]);

    int displacement = i * 8;
    
    checkCudaErrors(cudaMemcpyAsync(X_cord_slip + displacement   ,x, sizeof(double)*8,cudaMemcpyHostToDevice,FE3d_stream_slip));
    checkCudaErrors(cudaMemcpyAsync(Y_cord_slip + displacement   ,y, sizeof(double)*8,cudaMemcpyHostToDevice,FE3d_stream_slip));
    checkCudaErrors(cudaMemcpyAsync(Z_cord_slip + displacement   ,z, sizeof(double)*8,cudaMemcpyHostToDevice,FE3d_stream_slip));
    // cudaDeviceSynchronize();
    delete[] x; delete[] y; delete[] z;
  }
}

#endif  // #Ifdef _CUDA


////////////////////////// ARCHIVED CODES SECTION ///////////////////////////////////////////////////////////////
/*

void FE3D_ALE::pickDOFsOfFreeSlipBoundaries(TFESpace3D* gridfespace, std::vector<int> freeSlipBoundIds,std::vector<int> boundIds)
{
    int N_cells = coll->GetN_Cells();
    TBaseCell* currentCell;
    int MaxLen;
    int N_Joints;
    const int* TmpLen; const int* TmpFV;
    BoundCond Bdcond;
    TBoundFace* Bdface;					// Pointer to Boundary Face in 3D Cell
	TBoundComp *BoundComp;
    TVertex* currentVertex;
    bool cell_setflag = false;
    int* GlobalNumbers;
    int* BeginIndex;
    int N_Movfaces = 0;
    GlobalNumbers = gridfespace->GetGlobalNumbers();
	BeginIndex = gridfespace->GetBeginIndex();	

	cout << " Function - Pick Free SLIP BD's " <<endl;
    // Member Variables  -- Do  not Initalise them 
    N_bd_FreeSlip_Vertex = 0;
    N_bd_FreeSlip_Cells = 0;
    N_bd_FreeSlip_Joints = 0;
	N_bd_EdgeFreeSlip_Vertex = 0;
    N_bd_EdgeFreeSlip_Cells = 0;
    N_bd_EdgeFreeSlip_Joints = 0;


    for (int cellNr = 0 ; cellNr < N_cells ; cellNr++)
    {
        currentCell = coll->GetCell(cellNr);
        //get the Joints and the local numbering of Joints
        int* GlobalDOF = GlobalNumbers + BeginIndex[cellNr];
		FE3D elementId = gridfespace->GetFE3D(cellNr, currentCell);
        TFE3D *element = TFEDatabase3D::GetFE3D(elementId);
        TFEDesc3D *fedesc = element->GetFEDesc3D();
        currentCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        N_Joints = currentCell->GetN_Joints();
		int JointsPerCell = 0;
		currentCell->SetClipBoard(5);
        for ( int jointId = 0 ; jointId < N_Joints ; jointId++)  // Joints refer to Faces in 3D
		{
			TJoint* Joint = currentCell->GetJoint(jointId);
			if(Joint->GetType() == BoundaryFace)
			{
				// cout << "  Boundary Face " <<endl;
                Bdface = (TBoundFace*)Joint;
	  			BoundComp = Bdface->GetBoundComp();
	  			int bdid = BoundComp->GetID();
                std::vector<int>::iterator it = std::find(freeSlipBoundIds.begin(), freeSlipBoundIds.end(), bdid);
				if(it != freeSlipBoundIds.end())                        
				{
                    N_freeSurfaceJoints++;
					if(cell_setflag == FALSE){
						N_freeSurfaceCells++;					
					}   
					int *JointDOF = fedesc->GetJointDOF(jointId);
					int N_Vertices = TmpLen[jointId];
					N_bd_FreeSlip_Joints++;
					JointsPerCell++;
					cout << " Cell : " <<cellNr << " Joint : " <<jointId <<" BD ID : " <<bdid << " JPC : " <<JointsPerCell <<endl;
					//freeSurfaceJoints.emplace_back(jointId);
					for ( int vert = 0 ; vert < fedesc->GetN_JointDOF() ; vert++)
					{
						int local_vertex_no = TmpFV[jointId*MaxLen + vert];
						currentVertex = currentCell->GetVertex(vert) ;
						int glob_vertex_no = GlobalDOF[JointDOF[vert]];
						auto findIndex = std::find(Bd_FreeSlip_Vertex.begin(), Bd_FreeSlip_Vertex.end(), glob_vertex_no);
						// if (std::find(Bd_FreeSlip_Vertex.begin(), Bd_FreeSlip_Vertex.end(), glob_vertex_no) == Bd_FreeSlip_Vertex.end()) 
                        // {
                            if (std::find(freeSurfaceVertex.begin(), freeSurfaceVertex.end(), glob_vertex_no) == freeSurfaceVertex.end()) 
							{
								if(JointsPerCell == 1 && currentCell->GetClipBoard() == 5)  // Joint Comming for first Time  - FACE DOF
								{
									Bd_FreeSlip_Vertex.push_back(glob_vertex_no);
									Bd_FreeSlip_VertexLocal.emplace_back(local_vertex_no);
									Bd_FreeSlip_Cells.emplace_back(cellNr);
									Bd_FreeSlip_Joints.emplace_back(jointId);
									N_bd_FreeSlip_Vertex++;
									cell_setflag = TRUE ;
								}
								else      // Coming for Multiple Time - Edge DOF
								{
									Bd_EdgeFreeSlip_Vertex.push_back(glob_vertex_no);
									Bd_EdgeFreeSlip_VertexLocal.emplace_back(local_vertex_no);
									Bd_EdgeFreeSlip_Cells.emplace_back(cellNr);
									Bd_EdgeFreeSlip_Joints.emplace_back(jointId);
									N_bd_EdgeFreeSlip_Vertex++;
								}
								
                            }
						// }
					}
				}
			}
		
        }
    }

	for ( int j =  0 ; j < Bd_FreeSlip_Vertex.size() ; j++ )	
		cout <<  Bd_FreeSlip_Vertex[j] << " " ;
	cout<<endl;
	for ( int j =  0 ; j < Bd_EdgeFreeSlip_Vertex.size() ; j++ )	
		cout <<  Bd_EdgeFreeSlip_Vertex[j] << " " ;
		
    cout << " Number of FreeSlip Boundary cells** : " << N_bd_FreeSlip_Vertex <<endl;
	cout << " Number of FreeSlip Edge Boundary cells : " << N_bd_EdgeFreeSlip_Vertex <<endl;

} 

*/

/*
 
? Archived for Impose free Slip boundary as per Sir's code 


void FE3D_ALE::impose_FreeSlip_BoundaryCondition(TSquareStructure3D *sqstructure,double* a11, double* a12,double* a13,
												double* a21, double* a22, double* a23,
												double* a31, double* a32, double* a33 , double* rhs,int length)
{
	int* RowPtr = sqstructure->GetRowPtr();
	int* ColPtr = sqstructure->GetKCol();

	double Row1_val1 = 0,Row1_val2 = 0,Row1_val3 = 0,Row2_val1 = 0,
		   Row2_val2 = 0, Row2_val3 = 0,Row3_val1 = 0,Row3_val2 = 0,Row3_val3 = 0;

	for ( int vertexNo = 0 ; vertexNo < N_bd_FreeSlip_Vertex ; vertexNo++)
	{
		// Replace the 1st row with the Normal Values Directly
		int DOF = Bd_FreeSlip_Vertex[vertexNo];
		int Begin = RowPtr[DOF];
		int End   =  RowPtr[DOF + 1];
		
		for ( int rowIndex = Begin ; rowIndex < End ; rowIndex ++ )
		{
			// First Column of Blocks
			double v1 = a11[rowIndex];
      		double v2 = a21[rowIndex];
      		double v3 = a31[rowIndex];
			 cout << "[v1,v2,v3] : ["<<v1<<","<<v2<<","<<v3<<"]"<<endl;
			 cout << "[t1,t2,t3] : ["<<Bd_TangentA_1[vertexNo]<<","<<Bd_TangentA_2[vertexNo]<<","<<Bd_TangentA_3[vertexNo]<<"]"<<endl;
			a11[rowIndex] = Bd_TangentA_1[vertexNo]*v1 + Bd_TangentA_2[vertexNo]*v2 + Bd_TangentA_3[vertexNo]*v3;
			a21[rowIndex] = Bd_TangentB_1[vertexNo]*v1 + Bd_TangentB_2[vertexNo]*v2 + Bd_TangentB_3[vertexNo]*v3;
			cout << " a11,a21 : [" << a11[rowIndex] << "," << a21[rowIndex] <<"]" <<endl;
			// second column of blocks
			v1 = a12[rowIndex];
			v2 = a22[rowIndex];
			v3 = a32[rowIndex];
			
			a12[rowIndex] = Bd_TangentA_1[vertexNo]*v1 + Bd_TangentA_2[vertexNo]*v2 + Bd_TangentA_3[vertexNo]*v3;
			a22[rowIndex] = Bd_TangentB_1[vertexNo]*v1 + Bd_TangentB_2[vertexNo]*v2 + Bd_TangentB_3[vertexNo]*v3;

			
			// third column of blocks
			v1 = a13[rowIndex];
			v2 = a23[rowIndex];
			v3 = a33[rowIndex];
			
			a13[rowIndex] = Bd_TangentA_1[vertexNo]*v1 + Bd_TangentA_2[vertexNo]*v2 + Bd_TangentA_3[vertexNo]*v3;
			a23[rowIndex] = Bd_TangentB_1[vertexNo]*v1 + Bd_TangentB_2[vertexNo]*v2 + Bd_TangentB_3[vertexNo]*v3;
			

			if(ColPtr[rowIndex] == DOF)
			{
				a31[rowIndex] = Bd_normal_1[vertexNo];
        		a32[rowIndex] = Bd_normal_2[vertexNo];
        		a33[rowIndex] = Bd_normal_3[vertexNo];
			}

			else
			{
				a31[rowIndex] = 0;
        		a32[rowIndex] = 0;
        		a33[rowIndex] = 0;
			}	
		}
		double rhs1 = rhs[DOF];
		double rhs2 = rhs[DOF+length];
		double rhs3 = rhs[DOF + length];
		rhs[DOF] = Bd_TangentA_1[vertexNo]*rhs1 + Bd_TangentA_2[vertexNo]*rhs1 + Bd_TangentA_3[vertexNo]*rhs1;
		rhs[DOF+length] = Bd_TangentB_1[vertexNo]*rhs1 + Bd_TangentB_2[vertexNo]*rhs1 + Bd_TangentB_3[vertexNo]*rhs1;
		rhs[DOF + 2*length] = 0.0;
	}

	// FOR EDGE DOF 
	for ( int vertexNo = 0 ; vertexNo < N_bd_EdgeFreeSlip_Vertex ; vertexNo++)
	{
		// Replace the 1st row with the Normal Values Directly
		int DOF = Bd_EdgeFreeSlip_Vertex[vertexNo];
		int Begin = RowPtr[DOF];
		int End   =  RowPtr[DOF + 1];
		
		for ( int rowIndex = Begin ; rowIndex < End ; rowIndex ++ )
		{
			// First Column of Blocks
			double v1 = a11[rowIndex];
      		double v2 = a21[rowIndex];
      		double v3 = a31[rowIndex];
			 // cout << "[v1,v2,v3] : ["<<v1<<","<<v2<<","<<v3<<"]"<<endl;
			a11[rowIndex] = Bd_edge_TangentA_1[vertexNo]*v1 + Bd_edge_TangentA_2[vertexNo]*v2 + Bd_edge_TangentA_3[vertexNo]*v3;

			// cout << " a11,a21 : [" << a11[rowIndex] << "," << a21[rowIndex] <<"]" <<endl;
			// second column of blocks
			v1 = a12[rowIndex];
			v2 = a22[rowIndex];
			v3 = a32[rowIndex];
			
			a12[rowIndex] = Bd_edge_TangentA_1[vertexNo]*v1 + Bd_edge_TangentA_2[vertexNo]*v2 + Bd_edge_TangentA_3[vertexNo]*v3;
			//a22[rowIndex] = Bd_TangentB_1[vertexNo]*v1 + Bd_TangentB_2[vertexNo]*v2 + Bd_TangentB_3[vertexNo]*v3;

			
			// third column of blocks
			v1 = a13[rowIndex];
			v2 = a23[rowIndex];
			v3 = a33[rowIndex];
			
			a13[rowIndex] = Bd_edge_TangentA_1[vertexNo]*v1 + Bd_edge_TangentA_2[vertexNo]*v2 + Bd_edge_TangentA_3[vertexNo]*v3;
			// a23[rowIndex] = Bd_TangentB_1[vertexNo]*v1 + Bd_TangentB_2[vertexNo]*v2 + Bd_TangentB_3[vertexNo]*v3;
			

			if(ColPtr[rowIndex] == DOF)
			{
				a21[rowIndex] = Bd_edge_normalA_1[vertexNo];
        		a22[rowIndex] = Bd_edge_normalA_2[vertexNo];
        		a23[rowIndex] = Bd_edge_normalA_3[vertexNo];

				a31[rowIndex] = Bd_edge_normalB_1[vertexNo];
        		a32[rowIndex] = Bd_edge_normalB_2[vertexNo];
        		a33[rowIndex] = Bd_edge_normalB_3[vertexNo];
			}

			else
			{
				a21[rowIndex] = 0.;
        		a22[rowIndex] = 0.;
        		a23[rowIndex] = 0.;

				a31[rowIndex] = 0.;
        		a32[rowIndex] = 0.;
        		a33[rowIndex] = 0.;
			}	
		}
		double rhs1 = rhs[DOF];
		double rhs2 = rhs[DOF+length];
		double rhs3 = rhs[DOF + length];
		rhs[DOF] = Bd_edge_TangentA_1[vertexNo]*rhs1 + Bd_edge_TangentA_2[vertexNo]*rhs1 + Bd_edge_TangentA_3[vertexNo]*rhs1;
		rhs[DOF+length] = 0.0;
		rhs[DOF + 2*length] = 0.0;

	}


}


// Extracting free surface mesh velocity 
		// double* u_values = new double[3]();
		// double* meshVelocityValues = new double[3]();
		// //Get the Global Co-ordinates of the DOF in Velocity FESpace
		// TBaseCell* VelocityCell =  MovCells[i];
		// int freeSurface_globalID = freeSurfaceVertex[i];
		// switch(RefTransVelocity)
    	// {
		// 	case HexaTrilinear: 
		// 	{
		// 		int*  BeginIndex = Velocity_FEvect->GetFESpace3D()->GetBeginIndex();
		// 		int* GlobalNumbers = Velocity_FEvect->GetFESpace3D()->GetGlobalNumbers();

		// 		FE3D FE_ID = Velocity_FEvect->GetFESpace3D()->GetFE3D(freeSurfaceCells[vertex_number], currentCell);
		// 		TFE3D* FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
		// 		RefTrans = FE_Obj->GetRefTransID();
		// 		TFEDatabase3D::SetCellForRefTrans(currentCell, RefTrans);
		// 		// ------ Calculate the Values of the Velocity Function at the Nodal Points   ------- //
		// 		bf = ele_Velocity->GetBaseFunct3D();
		// 		int N_BaseFunct = bf->GetDimension();

		// 		double* uref = new double[N_BaseFunct];

		// 		// bf->GetDerivatives(D000, xiq, etaq, zetaq, uref);
		// 		bf->GetDerivatives(D000,xi_freeSurfaceNodes[vertex_number],
		// 			 eta_freeSurfaceNodes[vertex_number], zeta_freeSurfaceNodes[vertex_number], uref);


		// 		double u = 0,val;
		// 		int* Numbers = GlobalNumbers + BeginIndex[freeSurfaceCells[vertex_number]];
				
		// 		for ( int k = 0 ; k < 3 ; k++)
		// 		{
		// 			double* Values = Velocity_FEvect->GetComponent(k)->GetValues();
		// 			for(int j=0;j<N_BaseFunct;j++)
		// 			{
		// 				val = Values[Numbers[j]];
		// 				// cout << j << " " << val << endl;
		// 				meshVelocityValues[k]  +=  uref[j]*val;
		// 			} 
		// 		}
		// 		delete[] uref;
		// 		break;
		// 	}

		// 	case TetraAffin:
		// 	{
		// 		int*  BeginIndex = Velocity_FEvect->GetFESpace3D()->GetBeginIndex();
		// 		int* GlobalNumbers = Velocity_FEvect->GetFESpace3D()->GetGlobalNumbers();

		// 		FE3D FE_ID = Velocity_FEvect->GetFESpace3D()->GetFE3D(freeSurfaceCells[vertex_number], currentCell);
		// 		TFE3D* FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
		// 		RefTrans = FE_Obj->GetRefTransID();
		// 		TFEDatabase3D::SetCellForRefTrans(currentCell, RefTrans);
		// 		// ------ Calculate the Values of the Velocity Function at the Nodal Points   ------- //
		// 		bf = ele_Velocity->GetBaseFunct3D();
		// 		int N_BaseFunct = bf->GetDimension();

		// 		double* uref = new double[N_BaseFunct];

		// 		// bf->GetDerivatives(D000, xiq, etaq, zetaq, uref);
		// 		bf->GetDerivatives(D000,xi_freeSurfaceNodes[vertex_number],
		// 			 eta_freeSurfaceNodes[vertex_number], zeta_freeSurfaceNodes[vertex_number], uref);


		// 		double u = 0,val;
		// 		int* Numbers = GlobalNumbers + BeginIndex[freeSurfaceCells[vertex_number]];
				
		// 		for ( int k = 0 ; k < 3 ; k++)
		// 		{
		// 			double* Values = Velocity_FEvect->GetComponent(k)->GetValues();
		// 			for(int j=0;j<N_BaseFunct;j++)
		// 			{
		// 				val = Values[Numbers[j]];
		// 				// cout << j << " " << val << endl;
		// 				meshVelocityValues[k]  +=  uref[j]*val;
		// 			} 
		// 		}

		// 		delete[] uref;
		// 		break;
		// 	}
			


		// 	default:
		// 	{
		// 		cout << " ERROR TYPE : Method Not yet Implemented "<<endl;
		// 		cout << " ERROR DESC : Method to Calculate the Normal has not been Implemented for REFRTANS3D type : "<< referenceTransformation <<endl;
		// 		cout << " ERROR LOCATION : Class DeformMesh3D , Function: get_surface_normals " <<endl;
		// 		exit(0);
		// 	}
		// }

		// // Velocity_FEvect->FindValueLocal(MovCells[i], freeSurfaceCells[i], 
        //                             //   x_coord,y_c		// 
		// 							//   Mesh_FESpace->GetDOFPosition(freeSurfaceVertex[i], x_coord,y_coord, z_coord);
		// // double uVal[4] ={0}, vVal[4] = {0} , wVal[4] = {0};
		// // Velocity_FEvect->GetComponent(0)->FindGradientLocal(MovCells[i],freeSurfaceCells[i],x_coord, y_coord, z_coord, uVal);
		// // Velocity_FEvect->GetComponent(1)->FindGradientLocal(MovCells[i],freeSurfaceCells[i],x_coord, y_coord, z_coord, vVal);
		// // Velocity_FEvect->GetComponent(2)->FindGradientLocal(MovCells[i],freeSurfaceCells[i],x_coord, y_coord, z_coord, wVal);

		// // int checker = 0; 

		// // for ( int k = 0 ; k < 3 ; k++)
		// // {
		// // 	if ( (fabs( uVal[0] - meshVelocityValues[0]) > 1e-9 )
		// // 		|| (fabs( vVal[0] - meshVelocityValues[1]) > 1e-9 )
		// // 		|| (fabs( wVal[0] - meshVelocityValues[2]) > 1e-9 ))

		// // 	{
		// // 		checker++;
		// // 	}
		// // }

		// // if(checker) cout << " ERORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR" <<endl;oord,z_coord, 
        // //                               u_values);


				
		// u_values[0] = meshVelocityValues[0];
		// u_values[1] = meshVelocityValues[1];
		// u_values[2] = meshVelocityValues[2];


		// //////////// - END -- CODE BLOCK B - To get the Normal Vectors and Velocity Values at the nodal Points ///////////////////////
		// //cout << " Velocity Vector  : [ " << u_values[0] << ", " << u_values[1] << ", " << u_values[2] << " ] " <<endl;
		// double alpha = 0;
		// alpha = u_values[0]*norm_1 + u_values[1]*norm_2 + u_values[2]*norm_3;	
		// // cout << " Alpha : " << alpha <<endl;

		
	
		// // ASsign the ith component of Mesh Velocity as ( u1.n1 + u2.n2 + u3.n3 ) * Normal
        // // HARDCODED - For 
		// meshVelocityatFreeSurface_1[vertex_number] = alpha*norm_1;   // Value of w1 of Global DOF freeSurfaceVertex[i]
		// meshVelocityatFreeSurface_2[vertex_number] = alpha*norm_2;   // Value of w2 of Global DOF freeSurfaceVertex[i]
		// meshVelocityatFreeSurface_3[vertex_number] = alpha*norm_3;   // Value of w3 of Global DOF freeSurfaceVertex[i]
		
		// // meshVelocityatFreeSurface_1[vertex_number] = u_values[0];   // Value of w1 of Global DOF freeSurfaceVertex[i]
		// // meshVelocityatFreeSurface_2[vertex_number] = u_values[1];   // Value of w2 of Global DOF freeSurfaceVertex[i]
		// // meshVelocityatFreeSurface_3[vertex_number] = u_values[2];   // Value of w3 of Global DOF freeSurfaceVertex[i]


		// delete[] u_values;


*/