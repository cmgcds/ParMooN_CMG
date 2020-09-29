#ifndef __CG__
#define __CG__

#include <ItMethod.h>

/** iteration method */
class TCg : public TItMethod
{
  protected:

    /** arrays for cg depending on number of dof */
    double *r;
    double *z;
    double *p;
    double *Ap;
    //TODO eingefuegt fuer Polak-Ribiere
    double *r_last;
    //TODO eingefuegt fuer Polak-Ribiere
    //TODO eingefuegt fer flexible Paper
    double **d;
    //TODO eingefuegt fer flexible Paper

    /** matrices for flexible gmres depending on restart */

  public:
    /** constructor */
    TCg(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
      int n_aux, int N_Unknowns, int scalar);

    /** destructor */
    ~TCg();

    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol,
      double *rhs);
};
#endif
