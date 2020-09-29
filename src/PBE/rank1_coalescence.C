#include "rank1_coalescence.h"
#include <stdio.h>
#include <iostream>

int convolution(piecewise_linear f, piecewise_linear g, piecewise_linear w)
{
	extern double stepping; 

	f.create_hierarchical_structure();
	g.create_hierarchical_structure();
	w.create_hierarchical_structure();
	discont_function_p dfp_f = f.ToDiscontFuncP();
	discont_function_p dfp_g = g.ToDiscontFuncP();
	discont_function_p dfp_w = w.ToDiscontFuncP();

	for(int i = 0; i < f.L; i++)
	{
		discont_function_setup_points(dfp_f, i, f.grid[0], f.grid[0] + f.length[i], f.start_value + i, f.end_value + i);
		discont_function_setup_points(dfp_g, i, g.grid[0], g.grid[0] + g.length[i], g.start_value + i, g.end_value + i);
		discont_function_setup_points(dfp_w, i, w.grid[0], w.grid[0] + w.length[i], w.start_value + i, w.end_value + i);
	}
	func_p pf = convert_discont_function(dfp_f);
	func_p pg = convert_discont_function(dfp_g);
	discont_function_del(dfp_f);
	discont_function_del(dfp_g);

    func_p pw = func_new(f.L-1, 1);
	set_gridstructure(pw, dfp_w->intervals);
	func_p convolution_result = faltung_fepc(pf, pg, pw, stepping);

    func_del(pf);
    func_del(pg);
    func_del(pw);
      
    discont_function_p  tmp_dfp_w = convert_func(convolution_result, dfp_w->intervals);
	discont_function_del(dfp_w);
	dfp_w = tmp_dfp_w;

    	func_del(convolution_result);

	w.FromDiscontFuncP(dfp_w);

	return 0;
}

piecewise_linear_continuous pl_aggregation(piecewise_linear f, piecewise_linear g, piecewise_linear w_plus)
{
	convolution(f, g, w_plus);

	piecewise_linear_continuous w;
	w.IsContinuousProjectionOf(w_plus);

	int n = f.n, i;

	piecewise_linear_continuous fc;
	fc.IsContinuousProjectionOf(f);
	piecewise_linear_continuous gc;
	gc.IsContinuousProjectionOf(g);
	double Integral_of_g = gc.integral();

	for(i = 0; i < n; i++)
	{
		w.coeff[i] = 0.5 * w.coeff[i] - fc.coeff[i] * Integral_of_g;
	}

//	std::cout << "Mass of the result function before the projection is " << w.mass() << std::endl;
	w.IsMassConservedProjectionOf(w);
//	std::cout << "Mass of the result function after the projection is " << w.mass() << std::endl;

	return w;
}

