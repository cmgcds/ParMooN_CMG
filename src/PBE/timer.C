#include "timer.h"
#include <stdio.h>

void timer::start_timer()
{
	clk = clock();
	bEnabled = true;
	bPaused = false;
	t = 0;
}

void timer::pause_timer()
{
	if(!bEnabled)
		return;
	
	bPaused = true;
	t += (double)(clock() - clk) / CLOCKS_PER_SEC;
}

void timer::resume_timer()
{
	if(!bEnabled)
		return;
	bPaused = false;
	clk = clock();
}

void timer::stop_timer()
{
	bPaused = false;
	bEnabled = false;
	if(!bPaused)
		t += (double)(clock() - clk) / CLOCKS_PER_SEC;
}

double timer::get_time()
{
	if(!bEnabled)
		return t;
	return -1;
}