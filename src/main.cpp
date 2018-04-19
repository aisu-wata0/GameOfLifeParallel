
#include<iostream>
#include<unistd.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<sstream>

#include "Game.h"

using namespace G;

int main()
{
	Game g(40);
	
	g.init("Glider");
	
	usleep(1000000);
	
	for(long i = 0; true ; ++i)
	{
		g.print();
		g.nextGen();
		
		ostringstream flushStream;
		flushStream << "\r\x1b[" << g.size() << "A";
		//flushStream << "\r\x1bH";
		// \r - reset the poition to the beginning of the next line
		// \x1b - ESC
		// 	[numberA - move cursor up number lines
		// 	[numberB - move cursor down number lines
		// 	[numberC - move cursor right number chars
		// 	[numberD - move cursor left number chars
		// 	H - upper left corner
		cout << flushStream.str() << flush;
		usleep(75000);
	}

	return(0);
}