#include "complex.h"

idouble add(idouble a, idouble b) {
	idouble c; 
	c.m_dReal = a.m_dReal + b.m_dReal; 
	c.m_dImm = a.m_dImm + b.m_dImm; 
	return c; 
}

idouble subtract(idouble a, idouble b) {
	idouble c; 
	c.m_dReal = a.m_dReal - b.m_dReal; 
	c.m_dImm = a.m_dImm - b.m_dImm;
	return c; 
}

bool isZero(Complex p_dValue) {
	return (isZero(p_dValue.m_dReal) && isZero(p_dValue.m_dImm));
}
bool isReal(idouble a) {
	return isZero(a.m_dImm); 
}

idouble isqrt(Complex p_pIValue) {
	#define a p_pIValue.m_dReal
	#define b p_pIValue.m_dImm
	
	Complex res;
	
	double r = modulus(a, b);
	
	int eps;
	if (b != 0)
		eps = (int)sgn(b);
	else
		eps = 1;
	
	double alpha = (1/SQRT2);
	double dSqrtA = sqrt(r + a);
	double dSqrtB = sqrt(r - a);
	
	res.m_dReal = alpha * eps * dSqrtA;
	res.m_dImm = alpha * dSqrtB;
	
	return res;
	
	#undef a
	#undef b
}
