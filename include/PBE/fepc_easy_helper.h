
#ifndef _FEPC_EASY_HELPER
#define _FEPC_EASY_HELPER

#include "funktion.h"

func_p
func_multi(func_p f, func_p g);

func_p
func_div(func_p f, func_p g);

folgen_vektor_p
folgen_vektor_multi(folgen_vektor_p f, folgen_vektor_p g, int step);

folgen_vektor_p
folgen_vektor_div(folgen_vektor_p f, folgen_vektor_p g, int step);

fepc_real_t folge_norm2(folge_p folge, int step);

folge_p
folge_multi(folge_p f, folge_p g, int step);

folge_p
folge_div(folge_p f, folge_p g, int step);

vec_real_p get_mean_points(vec_p v, vec_p grad, int step);

#endif
