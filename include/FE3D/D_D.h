// mapper for hexahedron faces on both sides
static char D_D_Name[] = "D_D";
static char D_D_Desc[] = "discontinuous elements";
static int D_D_N0 = 0;
static int D_D_D = 0;
static int D_D_NPairs = 0;
static int *D_D_Pairs0 = NULL;
static int *D_D_Pairs1 = NULL;
static int *D_D_Pairs2 = NULL;
static int *D_D_Pairs3 = NULL;
static int *D_D_Pairs[4] = { (int *)D_D_Pairs0, (int *)D_D_Pairs1,
                               (int *)D_D_Pairs2, (int *)D_D_Pairs3 };

static int D_D_NNodes = 0;

static int D_D_NNoOpposite = 0;
static int **D_D_NoOpposite = NULL;

TFE3DMapper *D_D = new TFE3DMapper(D_D_Name, D_D_Desc,
                             D_D_N0, D_D_D,
                             D_D_NPairs, D_D_Pairs,
                             D_D_NNoOpposite, D_D_NoOpposite,
                             D_D_NNodes);
