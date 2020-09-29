#ifndef __BCGS__
#define __BCGS__

#include <ItMethod.h>

/** iteration method */
class TBcgs : public TItMethod
{
  protected:

  /** arrays for bcgs depending on number of dof */
	  double *r;
	  double *y;
	  double *p;
	  double *Ap;
	  double *v;
	  double *t;

  public:
    /** constructor */
    TBcgs(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
               int n_aux, int N_Unknowns, int scalar);

    /** destructor */
    ~TBcgs();

    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol,
                double *rhs);
};
#endif
