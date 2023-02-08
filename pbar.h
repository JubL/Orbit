#ifndef PBAR_H
#define PBAR_H

#include "obj.h"

#include <curses.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>
#include <sys/ioctl.h>
#include <stdio.h>

/**
	This class provides the progress bar feature.
*/
class Pbar
{
	public:
		Pbar();
		void pBarInit();
		void progressBar(double i);
		void showEstimatedTimeTillFinish(time_t& start, double advance, double settingsPosRow);
		void printInfo(double atmTemp, double soilTemp, Star& hoststar, Planet& exoplanet, std::string filename, time_t& start, double advance);
		void pBarStop();
	private:
		struct winsize _winsize;
		int32_t _lastCol;
};

#endif
