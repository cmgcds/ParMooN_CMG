#ifndef __FFT_FALTUNG_H
#define __FFT_FALTUNG_H

#include "config.h"

/*Berechnung der Faltung von zwei Arrays ueber die Fouriertransformation*/
fepc_real_t*
fft_faltung(fepc_real_t* a, int n_a, fepc_real_t* b, int n_b);

#endif
