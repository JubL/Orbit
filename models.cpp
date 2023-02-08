#include "models.h"

double temperature(double r_star, double distance, double t_star, double albedo, double epsilon)
{
	return t_star * sqrt(r_star / distance) * sqrt(sqrt((1 - albedo) / (4 * epsilon)));
}

double oneLayerTemperature(double r_star, double distance, double t_star, double albedo, double epsilon)
{
	return t_star * sqrt(r_star / distance) * sqrt(sqrt((1 - albedo) / (4 * (1 -  epsilon/2))));
}

double threeLayerTemperature(double r_star, double distance, double t_star, double albedo, double epsilon)
{
	return t_star * sqrt(r_star / distance) * sqrt(sqrt((1 - albedo) * ((1/pow((2 - epsilon), 2))+(epsilon/(4 * (2 - epsilon))))));
}

double cloudyTemperature(double r_star, double distance, double t_star, double albedo, double epsilon, double cloudRatio)
{
	return t_star * sqrt(r_star / distance) * sqrt(sqrt((1 - albedo)/(4 * epsilon * (1 - cloudRatio/2))));
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations" // this is neccessary to surpress compiler warnings of deprecated function calls inside this function (which itself is deprecated).
void manyLayers(std::vector<double>& temperatureArray,
		std::vector<double>& layerMass,
		std::vector<double>& numDens,
		Model& model,
		Star& hoststar,
		Planet& exoplanet,
		unsigned int it_beg,
		unsigned int it_end)
{
	if(it_end == 0) it_end = temperatureArray.size();
	if(it_beg > it_end) return;

	double crossSectionH2O = 2.22547e-25; // from exomol.com
	double crossSectionCO2 = 1.14906e-24;
	double crossSectionCH4 = 1.83999e-25;
	double crossSectionN2O = 0;
	double specHeatCapacitySurface = getSpecHeatCapacitySurface(model.getType());

	double power = 0;
	double starPower = 0;
	double deltaPower = 0;
	double absorbPower = 0;

	double numberDensityH2O = 0;
	double numberDensityCO2 = 0;
	double numberDensityCH4 = 0;
	double numberDensityN2O = 0;

	double temporaryTemp = 0;

	std::mutex mutex;

	for(unsigned int i = it_beg; i < it_end; i++) // i is the number of a layer
	{
		if(i == 0)
		{
			temporaryTemp = temperatureArray[0];
			power = radiateAndChangeTempRK4(temporaryTemp, 4 * M_PI * pow(exoplanet.getRadius() + 0.5 * model.getThicknessGasLayer(), 2), model.getDt(), layerMass[0], specHeatCapacitySurface);
			mutex.lock();
			temperatureArray[0] = temporaryTemp;
			mutex.unlock();

			starPower = model.getDt() * (1 - exoplanet.getAlbedo()) * pow(exoplanet.getRadius()/exoplanet.getDistance(), 2) * radiation(M_PI * pow(hoststar.getRadius(), 2), hoststar.getTemperature(), 1) * (1 - hoststar.getIRPortion());
			mutex.lock();
			temperatureArray[0] = newTemperature(starPower, specHeatCapacitySurface, layerMass[0], temperatureArray[0]);
			mutex.unlock();
		}
		else
		{
			temporaryTemp = temperatureArray[i];
			power = radiateAndChangeTempRK4(temporaryTemp, 4 * M_PI * pow(exoplanet.getRadius() + (i + 0.5) * model.getThicknessGasLayer(), 2), model.getDt(), layerMass[i], model.getSpecHeatCapacityAtm());
			mutex.lock();
			temperatureArray[i] = temporaryTemp;
			mutex.unlock();
		}

		deltaPower = power/2;
		for(unsigned int j = i + 1; j < temperatureArray.size(); j++) // for all layers above i
		{
			numberDensityH2O = numDens[j] * exoplanet.getVolumeRatioH2O();
			numberDensityCO2 = numDens[j] * exoplanet.getVolumeRatioCO2();
			numberDensityCH4 = numDens[j] * exoplanet.getVolumeRatioCH4();
			numberDensityN2O = numDens[j] * exoplanet.getVolumeRatioN2O();

			absorbPower = absorbtion(deltaPower, model.getThicknessGasLayer(), numberDensityH2O, crossSectionH2O, numberDensityCO2, crossSectionCO2, numberDensityCH4, crossSectionCH4, numberDensityN2O, crossSectionN2O);
			deltaPower -= absorbPower;
			mutex.lock();
			temperatureArray[j] = newTemperature(absorbPower, model.getSpecHeatCapacityAtm(), layerMass[j], temperatureArray[j]);
			mutex.unlock();
		}

		deltaPower = power/2;
		for(int k = i - 1; k >= 0; k--) // for all layers beneath i; this is not executed at i = 0 => k = -1
		{
			if(k == 0)
			{
				mutex.lock();
				temperatureArray[0] = newTemperature(deltaPower, specHeatCapacitySurface, layerMass[0], temperatureArray[0]); // if that layer is the surface absorb all radiation
				mutex.unlock();
			}
			else
			{
				numberDensityH2O = numDens[k] * exoplanet.getVolumeRatioH2O();
				numberDensityCO2 = numDens[k] * exoplanet.getVolumeRatioCO2();
				numberDensityCH4 = numDens[k] * exoplanet.getVolumeRatioCH4();
				numberDensityN2O = numDens[k] * exoplanet.getVolumeRatioN2O();

				absorbPower = absorbtion(deltaPower, model.getThicknessGasLayer(), numberDensityH2O, crossSectionH2O, numberDensityCO2, crossSectionCO2, numberDensityCH4, crossSectionCH4, numberDensityN2O, crossSectionN2O);
				deltaPower -= absorbPower;
				mutex.lock();
				temperatureArray[k] = newTemperature(absorbPower, model.getSpecHeatCapacityAtm(), layerMass[k], temperatureArray[k]);
				mutex.unlock();
			}
		}
	}
}
#pragma GCC diagnostic pop

void deepOcean(std::vector<double>& oceanTemp, std::vector<double>& oceanMass, std::vector<double>& temp, Planet& exoplanet, Model& model)
{
	oceanTemp[0] = temp[0];
	double power = 0;
	for(int i = 0; i < (int)oceanTemp.size() - 1; i++)
	{
		power = calcThermalConduction(oceanTemp, i, exoplanet, model); // returns the change in energy
		changeTemperature( power, getSpecHeatCapacitySurface(model.getType()), oceanMass[i],   oceanTemp[i]);
		changeTemperature(-power, getSpecHeatCapacitySurface(model.getType()), oceanMass[i+1], oceanTemp[i+1]);
	}
	temp[0] = oceanTemp[0];
}

void soil(std::vector<double>& soilTemp, std::vector<double>& soilMass, Star& hoststar, Planet& exoplanet, Model& model)
{
    starPower(model, hoststar, exoplanet, soilTemp[0], soilMass[0]);

	double power = 0;
	for(int i = 0; i < (int)soilTemp.size() - 1; i++)
	{
		power = calcThermalConduction(soilTemp, i, exoplanet, model); // returns the change in energy
		changeTemperature( power, getSpecHeatCapacitySurface(model.getType()), soilMass[i],   soilTemp[i]);
		changeTemperature(-power, getSpecHeatCapacitySurface(model.getType()), soilMass[i+1], soilTemp[i+1]);
	}
}

void twoStreamApproxMod(std::vector<double>& atmTemp, std::vector<double>& atmMass, std::vector<double>& numDens, Model& model, Star& hoststar, Planet& exoplanet)
{
	std::mutex mutex;

	double starPower = 0;
	double absorbPower = 0;

	double specHeatCapacitySurface = getSpecHeatCapacitySurface(model.getType());

	// This ist the energy of the IR portion of the stellar radiation
	starPower = model.getDt() * (1 - exoplanet.getAlbedo()) * pow(exoplanet.getRadius()/exoplanet.getDistance(), 2) * radiation(M_PI * pow(hoststar.getRadius(), 2), hoststar.getTemperature(), 1) * hoststar.getIRPortion();

	for(int i = atmTemp.size() - 1; i >= 0; i--)// let the IR radiation iterate through the atmosphere and deposit it's energy into the layers
	{
		if(i != 0)
		{
			absorbPower = absorb(starPower, model.getThicknessGasLayer(), numDens[i]);
			starPower -= absorbPower;
			changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[i], atmTemp[i]);
		}
		else
		{
			changeTemperature(starPower, specHeatCapacitySurface, atmMass[0], atmTemp[0]);
		}
	}
}

void oldGreenhouse(std::vector<double>& atmTemp,
                   std::vector<double>& atmMass,
                   std::vector<double>& numDens,
                   std::vector<double>& emitPower,
                   Model& model,
                   Star& hoststar,
                   Planet& exoplanet,
                   unsigned int it_beg,
                   unsigned int it_end)
{
	if(it_end == 0) it_end = atmTemp.size();
	if(it_beg > it_end) return;

	double specHeatCapacitySurface = getSpecHeatCapacitySurface(model.getType());

	double deltaPower = 0;
	double absorbPower = 0;

	if(it_beg == 0)
	{
		emitPower[0] = radiateAndChangeTempRK4(atmTemp[0], 4 * M_PI * pow(exoplanet.getRadius() + 0.5 * model.getThicknessGasLayer(), 2), model.getDt(), atmMass[0], specHeatCapacitySurface);

		starPower(model, hoststar, exoplanet, atmTemp[0], atmMass[0]);
	}
	else
	{
		emitPower[it_beg] = radiateAndChangeTempRK4(atmTemp[it_beg], 4 * M_PI * pow(exoplanet.getRadius() + (it_beg + 0.5) * model.getThicknessGasLayer(), 2), model.getDt(), atmMass[it_beg], model.getSpecHeatCapacityAtm());
	}

	for(unsigned int i = it_beg + 1; i < it_end; i++) // i is the number of a layer
	{
		emitPower[i] = radiateAndChangeTempRK4(atmTemp[i], 4 * M_PI * pow(exoplanet.getRadius() + (i + 0.5) * model.getThicknessGasLayer(), 2), model.getDt(), atmMass[i], model.getSpecHeatCapacityAtm());
	}

	for(unsigned int i = it_beg; i < it_end; i++) // i is the number of a layer
	{
		deltaPower = emitPower[i]/2;
		for(unsigned int j = i + 1; j < atmTemp.size(); j++) // for all layers above i
		{
			absorbPower = absorb(deltaPower, model.getThicknessGasLayer(), numDens[j]);
			deltaPower -= absorbPower;
			changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[j], atmTemp[j]);
		}

		deltaPower = emitPower[i]/2;
		for(int k = i - 1; k >= 0; k--) // for all layers beneath i; this is not executed at i = 0 => k = -1
		{
			if(k == 0)
			{
				changeTemperature(deltaPower, specHeatCapacitySurface, atmMass[0], atmTemp[0]); // if that layer is the surface absorb all radiation
			}
			else
			{
				absorbPower = absorb(deltaPower, model.getThicknessGasLayer(), numDens[k]);
				deltaPower -= absorbPower;
				changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[k], atmTemp[k]);
			}
		}
	}
}

void greenhouse(std::vector<double>& atmTemp,
                std::vector<double>& atmMass,
                std::vector<double>& soilTemp,
                std::vector<double>& soilMass,
                std::vector<double>& numDens,
                std::vector<double>& emitPower,
                Model& model,
                Planet& exoplanet)
{
	double specHeatCapacitySurface = getSpecHeatCapacitySurface(model.getType());

	double deltaPower = 0;
	double absorbPower = 0;

    // here the radiation from the surface is handled
    double soilPower = radiateAndChangeTempRK4(soilTemp[0], 4 * M_PI * pow(exoplanet.getRadius(), 2), model.getDt(), soilMass[0], specHeatCapacitySurface);
	for(unsigned int i = 0; i < atmTemp.size(); i++)
	{
		absorbPower = absorb(soilPower, model.getThicknessGasLayer(), numDens[i]);
		soilPower -= absorbPower;
		changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[i], atmTemp[i]);
	}

    // from here on the radiation from within the atmosphere is handled
	for(unsigned int i = 0; i < atmTemp.size(); i++) // all layers radiate power
	{
		emitPower[i] = radiateAndChangeTempRK4(atmTemp[i], 4 * M_PI * pow(exoplanet.getRadius() + (i + 0.5) * model.getThicknessGasLayer(), 2), model.getDt(), atmMass[i], model.getSpecHeatCapacityAtm());
	}

	for(unsigned int i = 0; i < atmTemp.size(); i++) // i is the number of a layer
	{
		deltaPower = emitPower[i]/2;    // TODO: Is this division by 2 really correct? Half the radiation goes upward, half goes downward, but we have calculated
		                                // the power for one surface only, and that makes 'emitPower' the power for that surface only. We need to calculate
		                                // the power for the other surface as well, instead of dividing by 2.
		for(unsigned int j = i + 1; j < atmTemp.size(); j++) // for all layers above i
		{
			absorbPower = absorb(deltaPower, model.getThicknessGasLayer(), numDens[j]);
			deltaPower -= absorbPower;
			changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[j], atmTemp[j]);
		}

		deltaPower = emitPower[i]/2;
		for(int k = i - 1; k >= -1; k--) // for all layers beneath i
		{
			if(k == -1)
			{
			    // if radiation hits the ground, absorb all of it
				changeTemperature(deltaPower*(1-exoplanet.getAlbedo()), specHeatCapacitySurface, soilMass[0], soilTemp[0]); // FIXME: ALBEDO?
				soilPower = deltaPower * exoplanet.getAlbedo();
            	for(unsigned int i = 0; i < atmTemp.size(); i++)
            	{
            		absorbPower = absorb(soilPower, model.getThicknessGasLayer(), numDens[i]);
            		soilPower -= absorbPower;
            		changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[i], atmTemp[i]);
            	}
			}
			else
			{
				absorbPower = absorb(deltaPower, model.getThicknessGasLayer(), numDens[k]);
				deltaPower -= absorbPower;
				changeTemperature(absorbPower, model.getSpecHeatCapacityAtm(), atmMass[k], atmTemp[k]);
			}
		}
	}
}
