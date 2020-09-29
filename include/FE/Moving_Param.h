#include <FEDatabase3D.h>

// ========================================================================
// parameter routine
// ========================================================================
void MovingParamsVelo3D(double *in, double *out);

// ========================================================================
// settings for parameter routine
// ========================================================================
int MovingN_FESpacesVelo = 3;
int MovingN_FctVelo = 6;
int MovingN_ParamFctVelo = 1;
int MovingN_FEValuesVelo = 6;
int MovingN_ParamsVelo = 3;
int MovingFEFctIndexVelo[6] = { 0, 1, 2, 3, 4, 5 };
MultiIndex3D MovingFEMultiIndexVelo[6] = { D000, D000, D000,
                                           D000, D000, D000 };
ParamFct *MovingFctVelo[1] = { MovingParamsVelo3D };
int MovingBeginParamVelo[1] = { 0 };

