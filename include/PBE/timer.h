#ifndef TIMER
#define TIMER

#include <time.h>

class timer
{
private:
	clock_t clk;
	bool bEnabled;
	bool bPaused;

	double t;

public:
	void start_timer();
	void pause_timer();
	void resume_timer();
	void stop_timer();

	double get_time();
};

#endif