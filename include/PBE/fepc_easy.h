#ifndef __FEPCEASY
#define __FEPCEASY

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "faltung.h"

typedef struct {
    fepc_real_t start;
    fepc_real_t end;
} interval_t; 

typedef interval_t * interval_p;

#ifndef INT_STEPS
#define INT_STEPS 10
#endif

#define SQRT_3 1.7320508075688772935274463415058723669428
#define DIV_1_6 0.16666666666666666666666666666666666666666666

typedef fepc_real_t (*Funcimpl_step) (int, int, int, fepc_real_t);

typedef fepc_real_t (*Funcimpl) (fepc_real_t);

/**
 * Sets the gridstructure in the needed form out of a common list of intervals.
 */
void set_gridstructure(func_p0_p function, interval_p* intervals, fepc_real_t stepping);
void set_gridstructure(func_p1_p function, interval_p* intervals, fepc_real_t stepping);

/**
 * Returns the h_l value corrensponding to the step.
 */
inline fepc_real_t get_h_l(int step, fepc_real_t& stepping) {
    return pow(.5, step)*stepping;
}

/**
 * Returns true if the vector is in the latter interval of the folge.
 */
inline bool_t is_in_latter_interval(int& v, folge_p0_p folge) {
	return v >= folge->lang / 2. ? FEPC_FALSE : FEPC_TRUE;
}

fepc_real_t get_value_at_step(func_p0_p function, fepc_real_t& x, int& step, fepc_real_t& stepping, int* count);

/**
 * Returns true if the vector is in the latter interval of the folge.
 */
inline bool_t is_in_latter_interval(int& v, folge_p1_p folge) {
	return v >= folge->lang / 2. ? FEPC_FALSE : FEPC_TRUE;
}

fepc_real_t get_value_at_step(func_p1_p function, fepc_real_t& x, int& step, fepc_real_t& stepping, int* count);

/**
 * Returns the value for phi_l (basis function) for the given vector, point x and step.
 */
fepc_real_t phi_l(int& step, int& v, int p, fepc_real_t& x, fepc_real_t& stepping);

#endif