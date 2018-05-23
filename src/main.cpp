
#include "main.hpp"

using namespace G;

int main(int argc, char **argv)
{
	omp_set_num_threads(4);
// Marker tutorial https://github.com/RRZE-HPC/likwid/wiki/TutorialMarkerC
#pragma omp parallel
{
	LIKWID_MARKER_THREADINIT;
}

	LIKWID_MARKER_INIT;

	srand(0);

	Game g(atoi(argv[1]));

	//g.init("Glider gun");
	g.init(0.25);

	for(long i = 0; i < 777 ; ++i)
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
