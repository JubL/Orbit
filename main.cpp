#include "constant.h"
#include "obj.h"
#include "models.h"
#include "functions.h"
#include "pbar.h"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <functional>
#include <chrono>
#include <typeinfo>
#include <signal.h>
#include <immintrin.h>
#include <csignal>

// TODO: in the manyLayers / greenhouse model: while looping throug the j and k loop: if the radiated energy is below x stop the loop
//       how much time is possibly saved? is this a problem considering the conservation of energy?
// TODO: improve the parallelisation, use OpenMP
// TODO: replace the vectors by pointers(?), introduce SIMD, align the data
// TODO: maybe use float instead of double?
// TODO: remove the -ncurse progress bar, replace it by a simple console output like in HCP/Serie_3/sinus.c
// TODO: work with energy instead of temperature to save a lot of calculations and thus time
// TODO: implement the surface-atmosphere coupling from Frierson et al. (2006) for heat transfer

// freezing of the water is disabled

int keyboardInterruptFlag = 0;

void signalHandler(int signum)
{
    std::cout << " Interrupt signal (" << signum << ") received." << std::endl;
    keyboardInterruptFlag = 1;
}

int calc(Model& model, Star& hoststar, Planet& exoplanet)
{
	std::string results;
//	results.reserve(7300000 * model.getCalcTime());

	std::vector<double> atmTemp;
	std::vector<double> atmMass;
	std::vector<double> oceanTemp;
	std::vector<double> oceanMass;
	std::vector<double> numDens;
	std::vector<double> emitPower;

	initVectors(model, hoststar, exoplanet, atmTemp, atmMass, oceanTemp, oceanMass, numDens, emitPower);

	long long unsigned int t = 0;
	double calculationTime = model.getCalcTime() * exoplanet.getOrbitalPeriod();
	int iterator = 0;

	exoplanet.setDistance();

	time_t start; // init and set the starting time
	time_t end;
	time(&start);

	bool isaTerminal = false;
	if(isatty(STDOUT_FILENO)) isaTerminal = true;
	Pbar progressBar; // initialize the progress bar
	if(isaTerminal)
	{
		initscr();
		progressBar.pBarInit();
		progressBar.progressBar(((double)t/calculationTime)*100);
		progressBar.printInfo(atmTemp[0], oceanTemp[0], hoststar, exoplanet, model.getFilename(), start, t/calculationTime);
	}

	while((double)t <= calculationTime)
	{
		exoplanet.setDistance(); // change these to planetMovement(t) ???
		exoplanet.calcPhaseAngle(model.getDt()); //
		t += model.getDt();
		iterator++;

//		oldGreenhouse(atmTemp, atmMass, numDens, emitPower, model, hoststar, exoplanet);
//		deepOcean(oceanTemp, oceanMass, atmTemp, exoplanet, model);
        greenhouse(atmTemp, atmMass, oceanTemp, oceanMass, numDens, emitPower, model, exoplanet);
        soil(oceanTemp, oceanMass, hoststar, exoplanet, model);
//		twoStreamApproxMod(atmTemp, atmMass, numDens, model, hoststar, exoplanet);
		// here we could introduce the method which executes atmospheric convection

		if(iterator == 100) // update the progress-bar and write data into the file every 100th timestep
		{
			if(keyboardInterruptFlag == 1) break;
			if(isaTerminal)
			{
				progressBar.progressBar(((double)t/calculationTime)*100);
				progressBar.printInfo(atmTemp[0], oceanTemp[0], hoststar, exoplanet, model.getFilename(), start, t/calculationTime);
			}
			results += std::to_string(t/model.getDt()) + " "
				+  std::to_string((t+model.getPreviousCalcTime()*exoplanet.getOrbitalPeriod())/86400.0) + " "
				+  std::to_string(exoplanet.getPhaseAngle()) + " "
				+  std::to_string(atmTemp[0]) + " "
				+  std::to_string(oceanTemp[0]) + " "
				+  std::to_string(exoplanet.getDistance() * std::cos(exoplanet.getPhaseAngle())) + " "
				+  std::to_string(exoplanet.getDistance() * std::sin(exoplanet.getPhaseAngle())) + "\n";
			iterator = 0;

			if(std::isnan(atmTemp[0]))
			{
				std::cout << "nan occured in atmTemp[0]" << std::endl;
				break;
			}
		}
	}
	if(isaTerminal)
	{
		progressBar.pBarStop();
	}
	time(&end);

    if(keyboardInterruptFlag == 1) std::cout << std::endl << "\033[1;31mTerminated by the user.\033[0m" << std::endl;

	double diff = difftime(end, start);
	writeDataToFile(model, hoststar, exoplanet, results, diff);

	std::cout << std::endl << "Required calculation time was " << diff << " seconds. With " << atmTemp.size() << " atmospheric layers and a time step of " << model.getDt() << " seconds." << std::endl;
	std::cout << "That are " << diff/60 << " minutes or " << diff/3600 << " hours." << std::endl;

	return 0;
}


int main()
{
    // register signal SIGINT and signal handler
    signal(SIGINT, signalHandler);

    #include "init.cpp"

/*	Sun hoststar;
	Earth exoplanet;

	exoplanet.setEccentricity(.4);

	Model model;

	model.setDt(200);
//	model.setCalcTime(10);
//	model.setThicknessSurfaceLayer(.05);
//	model.setOceanTempSize(200);
	model.setFilename("test.dat");
*/
	calc(model, hoststar, exoplanet);

	return 0;
}
