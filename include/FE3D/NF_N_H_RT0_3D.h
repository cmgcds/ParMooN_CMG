// ***********************************************************************
// Raviart-Thomas element of zero-th order on hexahedra, 3D
// ***********************************************************************

/* for all functionals */
static double NF_N_H_RT0_3D_Xi[]   = { 0, 0, 1, 0,-1, 0 };
static double NF_N_H_RT0_3D_Eta[]  = { 0,-1, 0, 1, 0, 0 };
static double NF_N_H_RT0_3D_Zeta[] = {-1, 0, 0, 0, 0, 1 };

/* face 0                               0 */
static double NF_N_H_RT0_3D_F0_Xi[]   = { 0 };
static double NF_N_H_RT0_3D_F0_Eta[]  = { 0 };
static double NF_N_H_RT0_3D_F0_Zeta[] = {-1 };

/* face 1                               1 */
static double NF_N_H_RT0_3D_F1_Xi[]   = { 0 };
static double NF_N_H_RT0_3D_F1_Eta[]  = {-1 };
static double NF_N_H_RT0_3D_F1_Zeta[] = { 0 };

/* face 2                               2 */
static double NF_N_H_RT0_3D_F2_Xi[]   = { 1 };
static double NF_N_H_RT0_3D_F2_Eta[]  = { 0 };
static double NF_N_H_RT0_3D_F2_Zeta[] = { 0 };

/* face 3                               3 */
static double NF_N_H_RT0_3D_F3_Xi[]   = { 0 };
static double NF_N_H_RT0_3D_F3_Eta[]  = { 1 };
static double NF_N_H_RT0_3D_F3_Zeta[] = { 0 };

/* face 4                               4 */
static double NF_N_H_RT0_3D_F4_Xi[]   = {-1 };
static double NF_N_H_RT0_3D_F4_Eta[]  = { 0 };
static double NF_N_H_RT0_3D_F4_Zeta[] = { 0 };

/* face 5                               5 */
static double NF_N_H_RT0_3D_F5_Xi[]   = { 0 };
static double NF_N_H_RT0_3D_F5_Eta[]  = { 0 };
static double NF_N_H_RT0_3D_F5_Zeta[] = { 1 };

static double *NF_N_H_RT0_3D_XiArray[6] = { 
                        NF_N_H_RT0_3D_F0_Xi,
                        NF_N_H_RT0_3D_F1_Xi,
                        NF_N_H_RT0_3D_F2_Xi,
                        NF_N_H_RT0_3D_F3_Xi,
                        NF_N_H_RT0_3D_F4_Xi,
                        NF_N_H_RT0_3D_F5_Xi };

static double *NF_N_H_RT0_3D_EtaArray[6] = { 
                        NF_N_H_RT0_3D_F0_Eta,
                        NF_N_H_RT0_3D_F1_Eta,
                        NF_N_H_RT0_3D_F2_Eta,
                        NF_N_H_RT0_3D_F3_Eta,
                        NF_N_H_RT0_3D_F4_Eta,
                        NF_N_H_RT0_3D_F5_Eta };

static double *NF_N_H_RT0_3D_ZetaArray[6] = { 
                        NF_N_H_RT0_3D_F0_Zeta,
                        NF_N_H_RT0_3D_F1_Zeta,
                        NF_N_H_RT0_3D_F2_Zeta,
                        NF_N_H_RT0_3D_F3_Zeta,
                        NF_N_H_RT0_3D_F4_Zeta,
                        NF_N_H_RT0_3D_F5_Zeta };

static double NF_N_H_RT0_3D_T[1] = {0};// ???
static double NF_N_H_RT0_3D_S[1] = {0};// ???

void NF_N_H_RT0_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  Functionals[0] = -PointValues[12] * 4.0;
  Functionals[1] = -PointValues[7]  * 4.0;
  Functionals[2] =  PointValues[2]  * 4.0;
  Functionals[3] =  PointValues[9]  * 4.0;
  Functionals[4] = -PointValues[4]  * 4.0;
  Functionals[5] =  PointValues[17] * 4.0;
}

void NF_N_H_RT0_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
                           double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  #ifdef __3D__
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 4, length == {4,4,4,4}
  Cell->GetVertex(faceVertex[4*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[4*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[4*face + 2])->GetCoords(x2,y2,z2);
  #endif
  // compute measure of this face
  s = sqrt( POW((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + POW((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + POW((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = PointValues[0]*s;
}

static int NF_N_H_RT0_3D_N_AllFunctionals = 6;
static int NF_N_H_RT0_3D_N_PointsAll = 6;
static int NF_N_H_RT0_3D_N_FaceFunctionals[] = { 1, 1, 1, 1, 1, 1 }; 
static int NF_N_H_RT0_3D_N_PointsFace[] = { 1, 1, 1, 1, 1, 1 };      

TNodalFunctional3D *NF_N_H_RT0_3D_Obj = new TNodalFunctional3D
        (NF_N_H_RT0_3D, NF_N_H_RT0_3D_N_AllFunctionals,
         NF_N_H_RT0_3D_N_FaceFunctionals, NF_N_H_RT0_3D_N_PointsAll,
         NF_N_H_RT0_3D_N_PointsFace,
         NF_N_H_RT0_3D_Xi, NF_N_H_RT0_3D_Eta, NF_N_H_RT0_3D_Zeta,
         NF_N_H_RT0_3D_XiArray, NF_N_H_RT0_3D_EtaArray,
         NF_N_H_RT0_3D_ZetaArray,
         NF_N_H_RT0_3D_T, NF_N_H_RT0_3D_S,
         NF_N_H_RT0_3D_EvalAll, NF_N_H_RT0_3D_EvalFace);
