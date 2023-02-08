#include "functions.h"
#include "obj.h"

#include <iostream>
#include <math.h>
#include <fstream>

/**
	This function is to calculate the cross section of gases which travel through the earth atmosphere and have a given atmospheric transmissivity value at the surface.
	@param initialValue This is the value at the beginning of the calculation.
	@param endValue This value is to be reached at the end of the calculation.
	@param exoplanet This method needs the Earth object to get the parameters like earth radius and earth mass from it.
	@param element Determine for which element you want to do the calculation, H2O, CO2 or CH4.
 */
void calcCrossSection(double initialValue, double endValue, Earth& exoplanet, Model& model, std::string element = "")
{
	double value = initialValue;
	double numberDensity = 0;
	double crossSection = pow(10, -25); // this is just a random guess

	double surfaceTemperature = 273 + 15;

	double margin = 0.00000001;
	double stepSize = 9.0;
	double direction = 0; // 1 = last value was too big; -1 = too small; 0 = yet unknown

	bool loopCondition = true;

	try
	{
		while(loopCondition)
		{
			value = initialValue;
			for(int i = 11000; i >= 0; i--) // i is in meters
			{
				numberDensity = atmNumberDensity(model, exoplanet, surfaceTemperature, i);
    				 if(element == "H2O") numberDensity *= 0.0003;
				else if(element == "CO2") numberDensity *= 3.45 * pow(10, -6);
				else if(element == "CH4") numberDensity *= 1.72 * pow(10, -8);
				value *= exp(-numberDensity * crossSection * 1);
			}
			if(value < endValue * (1 + margin) && value > endValue * (1 - margin)) loopCondition = false;
			else if(value >= endValue * (1 + margin))
			{
				if(direction == -1) stepSize /= 2;
				direction = 1;
				crossSection *= 1 + stepSize;
			}
			else if(value <= endValue * (1 - margin))
			{
				if(direction == 1) stepSize /= 2;
				direction = -1;
				crossSection /= 1 + stepSize;
			}
		}
	}
	catch(std::exception e)
	{
		std::cerr << e.what() << std::endl;
	}
	std::cout << "The cross section for " << element << " is " << crossSection << " m^2." << std::endl;
}

int main()
{
	Earth exoplanet;
	Model model;
	model.setThicknessGasLayer(1);

	calcCrossSection(1.0, .919253, exoplanet, model, "H2O"); // from NIST data (average of all data points used)
	calcCrossSection(1.0, .980869, exoplanet, model, "CO2");
	calcCrossSection(1.0, .984269, exoplanet, model, "CH4");

//	calcCrossSection(1.0, .91557, exoplanet, "H2O"); // from NIST data (2 to ~22.3 µm)
//	calcCrossSection(1.0, .980036, exoplanet, "CO2");
//	calcCrossSection(1.0, .984894, exoplanet, "CH4");

//	calcCrossSection(1.0, .976883, exoplanet, "H2O"); // from NIST data (8 to 12 µm)
//	calcCrossSection(1.0, .987205, exoplanet, "CO2");
//	calcCrossSection(1.0, .991038, exoplanet, "CH4");

        return 0;
}
