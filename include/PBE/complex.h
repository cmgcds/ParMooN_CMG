#ifndef MPIMIS_COMPLEX
#define MPIMIS_COMPLEX

#include "complextools.h"
#include "basictools.h"

typedef class Complex complex;
typedef complex* pcomplex;
typedef complex idouble;
typedef idouble* pidouble;

class Complex
{
public:
	Complex() { m_dReal = 0; m_dImm = 0; };
	~Complex() {};

	double m_dReal;
	double m_dImm;
};

idouble add(idouble a, idouble b);
idouble subtract(idouble a, idouble b);  

bool isZero(Complex p_dValue);
bool isReal(idouble a);

Complex isqrt(Complex p_pIValue);

#endif // MPIMIS_COMPLEX
