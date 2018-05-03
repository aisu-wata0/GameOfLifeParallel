
#include<iostream>
#include<unistd.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<sstream>
#include <omp.h>

#ifdef LIKWID_PERFMON
	#include <likwid.h>
#else
	#define LIKWID_MARKER_INIT
	#define LIKWID_MARKER_THREADINIT
	#define LIKWID_MARKER_SWITCH
	#define LIKWID_MARKER_REGISTER(regionTag)
	#define LIKWID_MARKER_START(regionTag)
	#define LIKWID_MARKER_STOP(regionTag)
	#define LIKWID_MARKER_CLOSE
	#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

#include "Game.h"

using namespace G;

int main()
{
LIKWID_MARKER_INIT;
	
	Game g(2048);
#ifdef LIKWID_PERFMON
	printf("Likiwid\n");
#endif
	g.init("Glider gun");
	
	for(long i = 0; i < 400 ; ++i)
	{
		/*g.print();*/
		g.nextGenB();
		/*
		ostringstream flushStream;
		flushStream << "\r\x1b[" << g.size() << "A";
		*/

		//flushStream << "\r\x1bH";
		// \r - reset the poition to the beginning of the next line
		// \x1b - ESC
		// 	[numberA - move cursor up number lines
		// 	[numberB - move cursor down number lines
		// 	[numberC - move cursor right number chars
		// 	[numberD - move cursor left number chars
		// 	H - upper left corner
		
		/*cout << flushStream.str() << flush;
		usleep(75000);*/
	}

LIKWID_MARKER_CLOSE;
	return(0);
}
