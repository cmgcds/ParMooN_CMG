#ifndef MPIMIS_TOOLS
#define MPIMIS_TOOLS

#include "basictools.h"

#include <math.h>
#include <stdlib.h>

inline bool isZero(double p_dValue);
inline double modulus(double p_dA, double p_dB);

inline bool isZero(double p_dValue) { // Tells if a number is zero considering a threshold 
	return (abs(p_dValue) < ZERO_THRESHOLD ? true : false);
}
inline double modulus(double p_dA, double p_dB) { // Gives the modulus of a complex number: p_dA + p_dB * i
	return sqrt(sqr(p_dA) + sqr(p_dB));
}

#endif // MPIMIS_TOOLS
