// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order on tetrahedra, 3D
// ***********************************************************************

//point
static double BDDF1T_tp = 0.166666666666667;

static double NF_N_T_BDDF1_3D_Xi[]  = {
    BDDF1T_tp, 0.5+BDDF1T_tp, BDDF1T_tp,   //face 0
    BDDF1T_tp, BDDF1T_tp, 0.5+BDDF1T_tp,   //face 1
    BDDF1T_tp, 1-2*BDDF1T_tp, BDDF1T_tp,   //face 2
    0,0,0                             //face 3
};
static double NF_N_T_BDDF1_3D_Eta[] = {
    BDDF1T_tp, BDDF1T_tp, 0.5+BDDF1T_tp,   //face 0
    0,0,0,                            //face 1
    1-2*BDDF1T_tp, BDDF1T_tp, BDDF1T_tp,   //face 2
    BDDF1T_tp, 0.5+BDDF1T_tp, BDDF1T_tp    //face 3
};
static double NF_N_T_BDDF1_3D_Zeta[]= {
    0,0,0,                            //face 0
    BDDF1T_tp, 0.5+BDDF1T_tp, BDDF1T_tp,   //face 1
    BDDF1T_tp, BDDF1T_tp, 1-2*BDDF1T_tp,   //face 2
    BDDF1T_tp, BDDF1T_tp, 0.5+BDDF1T_tp    //face 3
};

/* face 0                               0 */
static double NF_N_T_BDDF1_3D_F0_Xi[]   = { BDDF1T_tp, 0.5+BDDF1T_tp, BDDF1T_tp };
static double NF_N_T_BDDF1_3D_F0_Eta[]  = { BDDF1T_tp, BDDF1T_tp, 0.5+BDDF1T_tp };
static double NF_N_T_BDDF1_3D_F0_Zeta[] = { 0,0,0 };

/* face 1                               1 */
static double NF_N_T_BDDF1_3D_F1_Xi[]   = { BDDF1T_tp, BDDF1T_tp, 0.5+BDDF1T_tp };
static double NF_N_T_BDDF1_3D_F1_Eta[]  = { 0,0,0 };
static double NF_N_T_BDDF1_3D_F1_Zeta[] = { BDDF1T_tp, 0.5+BDDF1T_tp, BDDF1T_tp };

/* face 2                               2 */
static double NF_N_T_BDDF1_3D_F2_Xi[]   = { BDDF1T_tp, 1-2*BDDF1T_tp, BDDF1T_tp };
static double NF_N_T_BDDF1_3D_F2_Eta[]  = { 1-2*BDDF1T_tp, BDDF1T_tp, BDDF1T_tp };
static double NF_N_T_BDDF1_3D_F2_Zeta[] = { BDDF1T_tp, BDDF1T_tp, 1-2*BDDF1T_tp };


/* face 3                               3 */
static double NF_N_T_BDDF1_3D_F3_Xi[]   = { 0,0,0 };
static double NF_N_T_BDDF1_3D_F3_Eta[]  = { BDDF1T_tp, 0.5+BDDF1T_tp, BDDF1T_tp };
static double NF_N_T_BDDF1_3D_F3_Zeta[] = { BDDF1T_tp, BDDF1T_tp, 0.5+BDDF1T_tp };

static double *NF_N_T_BDDF1_3D_XiArray[4] = {
                        NF_N_T_BDDF1_3D_F0_Xi,
                        NF_N_T_BDDF1_3D_F1_Xi,
                        NF_N_T_BDDF1_3D_F2_Xi,
                        NF_N_T_BDDF1_3D_F3_Xi };

static double *NF_N_T_BDDF1_3D_EtaArray[4] = {
                        NF_N_T_BDDF1_3D_F0_Eta,
                        NF_N_T_BDDF1_3D_F1_Eta,
                        NF_N_T_BDDF1_3D_F2_Eta,
                        NF_N_T_BDDF1_3D_F3_Eta };

static double *NF_N_T_BDDF1_3D_ZetaArray[4] = {
                        NF_N_T_BDDF1_3D_F0_Zeta,
                        NF_N_T_BDDF1_3D_F1_Zeta,
                        NF_N_T_BDDF1_3D_F2_Zeta,
                        NF_N_T_BDDF1_3D_F3_Zeta };

static double NF_N_T_BDDF1_3D_T[] = {-100};// ???
static double NF_N_T_BDDF1_3D_S[] = {-100};// ???

void NF_N_T_BDDF1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  //face 0
  Functionals[0]  = -PointValues[24];
  Functionals[1]  = -PointValues[25];
  Functionals[2]  = -PointValues[26];
  //face 1
  Functionals[3]  = -PointValues[15];
  Functionals[4]  = -PointValues[16];
  Functionals[5]  = -PointValues[17];
  //face 2
  Functionals[6]  =  PointValues[6]+PointValues[18]+PointValues[30];
  Functionals[7]  =  PointValues[7]+PointValues[19]+PointValues[31];
  Functionals[8]  =  PointValues[8]+PointValues[20]+PointValues[32];
  //face 3
  Functionals[9]  = -PointValues[9];
  Functionals[10] = -PointValues[10];
  Functionals[11] = -PointValues[11];
}

void NF_N_T_BDDF1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int face,
                           double *PointValues, double *Functionals)
{
  double s; // size of face
  double x0,x1,x2,y0,y1,y2,z0,z1,z2;
  #ifdef __3D__
  // find vertices of this face, then their coordinates
  const int *faceVertex, *length;
  int MaxLen;
  Cell->GetShapeDesc()->GetFaceVertex(faceVertex, length, MaxLen);
  // now MaxLen == 3, length == {3,3,3,3}
  Cell->GetVertex(faceVertex[3*face    ])->GetCoords(x0,y0,z0);
  Cell->GetVertex(faceVertex[3*face + 1])->GetCoords(x1,y1,z1);
  Cell->GetVertex(faceVertex[3*face + 2])->GetCoords(x2,y2,z2);
  #endif
  // compute measure of this face
  s = sqrt( POW((y1-y0)*(z2-z0) - (z1-z0)*(y2-y0),2)
          + POW((z1-z0)*(x2-x0) - (x1-x0)*(z2-z0),2)
          + POW((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0),2) );
  Functionals[0] = PointValues[0]*s;
  Functionals[1] = PointValues[1]*s;
  Functionals[2] = PointValues[2]*s;
}

static int NF_N_T_BDDF1_3D_N_AllFunctionals = 12;
static int NF_N_T_BDDF1_3D_N_PointsAll = 12;
static int NF_N_T_BDDF1_3D_N_FaceFunctionals[] = { 3, 3, 3, 3 };
static int NF_N_T_BDDF1_3D_N_PointsFace[] = { 3, 3, 3, 3 };

TNodalFunctional3D *NF_N_T_BDDF1_3D_Obj = new TNodalFunctional3D
        (NF_N_T_BDDF1_3D, NF_N_T_BDDF1_3D_N_AllFunctionals,
         NF_N_T_BDDF1_3D_N_FaceFunctionals, NF_N_T_BDDF1_3D_N_PointsAll,
         NF_N_T_BDDF1_3D_N_PointsFace,
         NF_N_T_BDDF1_3D_Xi, NF_N_T_BDDF1_3D_Eta, NF_N_T_BDDF1_3D_Zeta,
         NF_N_T_BDDF1_3D_XiArray, NF_N_T_BDDF1_3D_EtaArray,
         NF_N_T_BDDF1_3D_ZetaArray,
         NF_N_T_BDDF1_3D_T, NF_N_T_BDDF1_3D_S,
         NF_N_T_BDDF1_3D_EvalAll, NF_N_T_BDDF1_3D_EvalFace);
