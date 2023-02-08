#include "functions.h"

void calcMean(double& meanValue, double currentValue, int n)
{
	meanValue += currentValue / n;
}

double distance(double e, double a, double phi)
{
	return a * (1 - std::pow(e, 2))/(1 + e * std::cos(phi));
}

double nextAngle(double e, double a, double r, double p, double dt)
{
	return pow(a/r, 2) * std::sqrt(1 - std::pow(e, 2)) * (2 * M_PI / p) * dt;
}

double wiensDisplacementLaw(double temp)
{
	return 2.89776829 * std::pow(10,6) / temp;
}

double radiation(double area, double temperature, double epsilon)
{
	return Constant::stefanBoltzmann * std::pow(temperature, 4) * area * epsilon;
}

double radiateAndChangeTemp(double& temperature, double area, double dt, double layerMass, double specHeatCapacity)
{
	double power = 0;
	double actualPower = 0;
	const double stepSize = .1;
	double actualTemperature = temperature;
	for(double i = 0; i < dt; i += stepSize)
	{
		actualPower = stepSize * radiation(area, actualTemperature);
		changeTemperature(-actualPower, specHeatCapacity, layerMass, actualTemperature);
		if(actualTemperature < 0)
		{
			actualTemperature = 0;
			break;
		}
		power += actualPower;
	}
	temperature = actualTemperature;
	return power;
}

double radiateAndChangeTempRK4(double& temperature, double area, double dt, double layerMass, double specHeatCapacity)
{
	std::mutex mutex;

	const double sbConst = -(Constant::stefanBoltzmann * area) / (specHeatCapacity * layerMass);
	double formerTemperature = temperature;

	double k1 = sbConst * pow(temperature, 4);
	double k2 = sbConst * pow(temperature + .5 * dt * k1, 4);
	double k3 = sbConst * pow(temperature + .5 * dt * k2, 4);
	double k4 = sbConst * pow(temperature +      dt * k3, 4);

	mutex.lock();
	temperature += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
	mutex.unlock();

	return (formerTemperature - temperature) * specHeatCapacity * layerMass;
}

double newTemperature(double heat, double specHeatCapacity, double mass, double oldTemperature)
{
	return heat/(specHeatCapacity * mass) + oldTemperature;
}

void changeTemperature(double heat, double specHeatCapacity, double mass, double& temperature)
{
	std::mutex mutex;
	mutex.lock();
	temperature += heat/(specHeatCapacity * mass);
	mutex.unlock();
}

double absorbtion(double initialQuantity, double thickness, double numberDensity1, double crossSection1, double numberDensity2, double crossSection2, double numberDensity3, double crossSection3, double numberDensity4, double crossSection4)
{
	return -initialQuantity * expm1(- thickness * (numberDensity1 * crossSection1 + numberDensity2 * crossSection2 + numberDensity3 * crossSection3 + numberDensity4 * crossSection4));
}

double absorb(double initialQuantity, double thickness, double numDens)
{
	return -initialQuantity * expm1(- thickness * numDens);
}

double atmNumberDensity(Model& model, Planet& exoplanet, double surfaceTemperature, int layerNumber)
{
	double r = ((double)layerNumber - 1/2) * model.getThicknessGasLayer();
	r += exoplanet.getRadius();

	double base = (r * surfaceTemperature) / (exoplanet.getRadius() * (surfaceTemperature + model.getTemperatureGradient() * exoplanet.getRadius() - model.getTemperatureGradient() * r));
	double exponent = - (Constant::gravitationalConstant * exoplanet.getMass() * model.getAtmMeanMolarMass() * model.getTemperatureGradient()) / (Constant::gasConstant * std::pow(surfaceTemperature + model.getTemperatureGradient() * exoplanet.getRadius(), 2));
	double argOfExp = - (Constant::gravitationalConstant * exoplanet.getMass() * model.getAtmMeanMolarMass()) / (Constant::gasConstant * (surfaceTemperature + model.getTemperatureGradient() * exoplanet.getRadius())) * (r - exoplanet.getRadius())/(r * exoplanet.getRadius());

	return exoplanet.getNumberDensityAtSurface() * std::pow(base, exponent) * std::exp(argOfExp);
}

double surfaceAtmNumberDensity(double surfaceAtmPressure, double surfaceAtmTemperature)
{
	return surfaceAtmPressure / (Constant::boltzmannConstant * surfaceAtmTemperature);
}

double meanMolarMass(double volumeRatio1, double molarMass1, double volumeRatio2, double molarMass2, double volumeRatio3, double molarMass3, double volumeRatio4, double molarMass4, double volumeRatio5, double molarMass5)
{
	return volumeRatio1 * molarMass1 + volumeRatio2 * molarMass2 + volumeRatio3 * molarMass3 + volumeRatio4 * molarMass4 + volumeRatio5 * molarMass5;
}

double massAtmLayer(double molarMass, double numberDensity, double layerNumber, double r_planet, double thickness)
{
	return molarMass * numberDensity / Constant::avogadroConstant * 4 * M_PI * pow(r_planet + (layerNumber + 1/2) * thickness, 2) * thickness;
}

void createTempVector(std::vector<double>& temp, Star& hoststar, Planet& exoplanet, Model& model)
{
	temp.reserve(model.getAtmTempSize());
	// The subtrahend improves the calcuation in such a manner that it prevents an overshooting of the temperature curve at the beginning of the simulation.
	double h = hoststar.getTemperature() * sqrt(hoststar.getRadius() / exoplanet.getSemiMajorAxis()) * sqrt(sqrt((1 - exoplanet.getAlbedo()) / 4)) - 30;
	temp.push_back(h);
	for(int i = 1; i < model.getAtmTempSize(); i++)
	{
		h -= 0.0065;
		temp.push_back(h);
	}
}

void createMassVector(std::vector<double>& mass, std::vector<double>& temp, std::vector<double>& numDens, Planet& exoplanet, Model& model)
{
	mass.reserve(temp.size());
	mass.push_back(modelSurfaceMass(exoplanet, model));
	for(int i = 1; i < (int)temp.size(); i++)
	{
		mass.push_back(massAtmLayer(model.getAtmMeanMolarMass(), numDens[i], i, exoplanet.getRadius(), model.getThicknessGasLayer()));
	}
}

void createNumDensVector(std::vector<double>& numDens, std::vector<double>& temp, Planet& exoplanet, Model& model)
{
	numDens.reserve(temp.size());
	for(int i = 0; i < (int)temp.size(); i++)
	{
		if(i == 0) numDens.push_back(exoplanet.getNumberDensityAtSurface());
		else numDens.push_back(atmNumberDensity(model, exoplanet, temp[0], i));
	}
}

double modelSurfaceMass(Planet& exoplanet, Model& model, int layer)
{
	double volume = 4 * M_PI * pow(exoplanet.getRadius() - layer * model.getThicknessSurfaceLayer(), 2) * model.getThicknessSurfaceLayer();
	double density = 0;
	std::string type = model.getType();

         if(type == "water")     density = 1000;
	else if(type == "swamp")     density = 1500; // Source: ?
	else if(type == "earthlike") density = .7 * 1000 + .3 * 1141; // the last number is highly imprecise!
	else if(type == "rock")      density = 2650; // feldspars
	else if(type == "granite")   density = 2700; // Source: ?
	else if(type == "desert")    density = 1800; // sand and gravel (dry)
	else if(type == "ice")       density = 916.7;

	return volume * density;
}

double getSpecHeatCapacitySurface(std::string type)
{
	double temp = 0;

         if(type == "water")     temp = 4186; // Source: CRC Handbook of Chemistry and Physics
	else if(type == "swamp")     temp = 1; // ?
	else if(type == "earthlike") temp = 1; // ?
	else if(type == "rock")      temp = 1; // ? // feldspars
	else if(type == "granite")   temp = 790; // Source: https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html
	else if(type == "desert")    temp = 830; // sand, quartz; Source: https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html
	else if(type == "ice")       temp = 2093; // Source: https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html
	else if(type == "snow")      temp = 2106;

	return temp;
}

double calcThermalConduction(std::vector<double>& temp, int i, Planet& exoplanet, Model& model)
{
	if(i >= (int)(temp.size() - 1)) return 0;

	double thermalConductivity = 0;
         if(model.getType() == "water")   thermalConductivity = 0.592;
	else if(model.getType() == "granite") thermalConductivity = 2.8;
	else if(model.getType() == "desert")  thermalConductivity = 0.58; // sand
	else if(model.getType() == "ice")     thermalConductivity = 2.33;
	else if(model.getType() == "snow")    thermalConductivity = 0.16;

//	if(model.getType() == "water" && temp[i] < 273.15 && temp[i+1] < 273.15) thermalConductivity = 2.33; // tread as ice

	double area = 4 * M_PI * pow(exoplanet.getRadius() - (i+1) * model.getThicknessSurfaceLayer(), 2);
	return - thermalConductivity * area * model.getDt() * (temp[i] - temp[i+1]) / model.getThicknessSurfaceLayer();
}

void writeMetadata(Model& model, Star& hoststar, Planet& exoplanet, double diff, std::string message)
{
	std::fstream f(model.getFilename(), std::ios::out);
	f << "# Star properties: " << "Mass: " << hoststar.getMass(1) << " solar masses, " << "Radius: " << hoststar.getRadius(1) << " solar radii, " << "SurfaceTemperature: " << hoststar.getTemperature() << " Kelvin" << std::endl;
	f << "# Exoplanet properties: " << "Mass: " << exoplanet.getMass(1) << " earth masses, " << "Radius: " << exoplanet.getRadius(1) << " earth radii, " << "Albedo: " << exoplanet.getAlbedo() << ", " << "Orbital eccentricity: " << exoplanet.getEccentricity() << ", " << "Semi Major Axis: " << exoplanet.getSemiMajorAxis(1) << " AU, " << "Orbital period: " << exoplanet.getOrbitalPeriod(1) << " earth days, " << "Volumeratio H2O: " << exoplanet.getVolumeRatioH2O() << ", " << "VolumeratioCO2: " << exoplanet.getVolumeRatioCO2() << ", " << "VolumeratioCH4: " << exoplanet.getVolumeRatioCH4() << ", " << "VolumeratioN2O: " << exoplanet.getVolumeRatioN2O() << std::endl;
	f << "# Other properties: " << "Timestep: " << model.getDt() << " seconds, " << "Thickness of the Surface Layer: " << model.getThicknessSurfaceLayer() << " meter, " << "Thickness of the gas layers: " << model.getThicknessGasLayer() << " meter, " << "Mean molar mass of the atmospheric gas: " << model.getAtmMeanMolarMass() << " kilogram, " << "Surface type: " << model.getType() << ", " << "Temperature Gradient: " << model.getTemperatureGradient() << " Kelvin per meter, " << " Number of ocean layers: " << model.getOceanTempSize() << ", " << "Starting Temperature of the ocean: " << model.getOceanStartingTemp() << " Kelvin, " << "Calculation Time: " << model.getCalcTime() << " earth years, " << "Specific heat capacity of the atmospheric gas: " << model.getSpecHeatCapacityAtm() << " Joule per kilogram per Kelvin" << std::endl;
	f << "# The calculation took " << std::to_string(diff) << " seconds." << std::endl;
	f << message << std::endl;
	f << std::endl;
	f.close();
}

void noteStarFlux(double& starPower, double& albedo, double& planetRadius)
{
	std::fstream file("fluxData.txt", std::ios::app);
	file << starPower /((1-albedo) * M_PI * pow(planetRadius, 2)) << std::endl;
	file.close();
}

void writeContinuationData(std::vector<double>& temp, std::vector<double>& layerMass, std::vector<double>& oceanTemp, std::vector<double>& oceanMass, Model& model, Star& hoststar, Planet& exoplanet)
{
	std::fstream f1("cont/continuationData_temp.dat", std::ios::out);
	f1 << "# temp data" << std::endl;
	for(int i = 0; i < (int)temp.size(); i++)
	{
		f1 << temp[i] << std::endl;
	}
	f1.close();

	std::fstream f2("cont/continuationData_layerMass.dat", std::ios::out);
	f1 << "# layerMass data" << std::endl;
	for(int i = 0; i < (int)layerMass.size(); i++)
	{
		f2 << layerMass[i] << std::endl;
	}
	f2.close();

	std::fstream f3("cont/continuationData_oceanTemp.dat", std::ios::out);
	f1 << "# ocenaTemp data" << std::endl;
	for(int i = 0; i < (int)oceanTemp.size(); i++)
	{
		f3 << oceanTemp[i] << std::endl;
	}
	f3.close();

	std::fstream f4("cont/continuationData_oceanMass.dat", std::ios::out);
	f1 << "# oceanMass data" << std::endl;
	for(int i = 0; i < (int)oceanMass.size(); i++)
	{
		f4 << oceanMass[i] << std::endl;
	}
	f4.close();

	std::fstream f5("cont/continuationData_model.dat", std::ios::out);
	f5 << "# model parameters" << std::endl;
	f5 << model.getDt() << " This is the time step." << std::endl;
	f5 << model.getThicknessSurfaceLayer() << " This is the thickness of the surface layer." << std::endl;
	f5 << model.getThicknessGasLayer() << " This is the thickness of the gas layers." << std::endl;
	f5 << model.getAtmMeanMolarMass() << " This is the mean molar mass of the 'air' molecules." << std::endl;
	f5 << model.getTemperatureGradient() << " This temperature gradient describes the drop in temperature along the height." << std::endl;
//	f5 << model.getOceanTempSize() << " This is the number of ocean layers, including the surface layer." << std::endl;
//	f5 << model.getOceanStartingTemp() << " This is the temperatur of the ocean at the begin of the calculations" << std::endl;
	f5 << model.getSpecHeatCapacityAtm() << " This is the specific heat capacity of the atmosphere." << std::endl;
	f5 << model.getCalcTime() << " This is the time the calculations covered in earth years." << std::endl;
//	f5 << model.getContinuation() << " ########" << std::endl;
	f5 << model.getType() << " This string describes the material the surface consits of, e.g. 'water'." << std::endl;
	f5 << model.getFilename() << " This is the filename the data is saved into." << std::endl;
	f5.close();

	std::fstream f6("cont/continuationData_hoststar.dat", std::ios::out);
	f6 << "# star parameters" << std::endl;
	f6 << hoststar.getMass(1) << " This is the stars mass in solar masses." << std::endl;
	f6 << hoststar.getRadius(1) << " This is the stars radius in solar radii." << std::endl;
	f6 << hoststar.getTemperature() << " This is the stars effective temperature in kelvin." << std::endl;
	f6.close();

	std::fstream f7("cont/continuationData_exoplanet.dat", std::ios::out);
	f7 << "# planet parameters" << std::endl;
	f7 << exoplanet.getMass(1) << " This is the planets mass in earth masses." << std::endl;
	f7 << exoplanet.getRadius(1) << " This is the planets radius in earth radii." << std::endl;
	f7 << exoplanet.getAlbedo() << " This is the planets albedo." << std::endl;
	f7 << exoplanet.getEccentricity() << " This is the planets orbital eccentricity." << std::endl;
	f7 << exoplanet.getSemiMajorAxis(1) << " This is the planets orbital semi major axis in astronomical units." << std::endl;
	f7 << exoplanet.getDistance() << " This is the planets distance to it's hoststar in meters." << std::endl;
	f7 << exoplanet.getOrbitalPeriod(1) << " This is the planets orbital period in earth days." << std::endl;
	f7 << exoplanet.getPhaseAngle() << " This is the current phase angle the planet is at." << std::endl;
	f7 << exoplanet.getNumberDensityAtSurface() << " This is the planet atmospheres number density at the surface." << std::endl;
	f7 << exoplanet.getVolumeRatioH2O() << " This is the planet atmospheres volume ration of H2O." << std::endl;
	f7 << exoplanet.getVolumeRatioCO2() << " This is the planet atmospheres volume ration of CO2." << std::endl;
	f7 << exoplanet.getVolumeRatioCH4() << " This is the planet atmospheres volume ration of CH4." << std::endl;
	f7 << exoplanet.getVolumeRatioN2O() << " This is the planet atmospheres volume ration of N2O." << std::endl;
	f7.close();
}

int loadContinuationData(std::vector<double>& temp, std::vector<double>& layerMass, std::vector<double>& oceanTemp, std::vector<double>& oceanMass, Model& model, Star& hoststar, Planet& exoplanet)
{
	std::ifstream f1("cont/continuationData_temp.dat");
	std::ifstream f2("cont/continuationData_layerMass.dat");
	std::ifstream f3("cont/continuationData_oceanTemp.dat");
	std::ifstream f4("cont/continuationData_oceanMass.dat");
	std::ifstream f5("cont/continuationData_model.dat");
	std::ifstream f6("cont/continuationData_hoststar.dat");
	std::ifstream f7("cont/continuationData_exoplanet.dat");

	if(!f1.good() || !f2.good() || !f3.good() || ! f4.good() || !f5.good() || !f6.good() || !f7.good()) return -1; // could be more precise

	std::string tempBuffer;

	try
	{
		int num = numberOfLines(f1) - 1;
		temp.reserve(num);
		while(std::getline(f1, tempBuffer))
		{
			if(tempBuffer[0] == '#') continue;
			temp.push_back(std::stod(tempBuffer));
		}
		f1.close();
	}
	catch (const std::exception& e)
	{
		f1.close();
		std::cout << "Initialisation of the 'temp' vector failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	try
	{
		int num = numberOfLines(f2) - 1;
		layerMass.reserve(num);
		while(std::getline(f2, tempBuffer))
		{
			if(tempBuffer[0] == '#') continue;
			layerMass.push_back(std::stod(tempBuffer));
		}
		f2.close();
	}
	catch (const std::exception& e)
	{
		f2.close();
		std::cout << "Initialisation of the 'layerMass' vector failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	try
	{
		int num = numberOfLines(f3) - 1;
		oceanTemp.reserve(num);
		while(std::getline(f3, tempBuffer))
		{
			if(tempBuffer[0] == '#') continue;
			oceanTemp.push_back(std::stod(tempBuffer));
		}
		f3.close();
	}
	catch (const std::exception& e)
	{
		f3.close();
		std::cout << "Initialisation of the 'oceanTemp' vector failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	try
	{
		int num = numberOfLines(f4) - 1;
		oceanMass.reserve(num);
		while(std::getline(f4, tempBuffer))
		{
			if(tempBuffer[0] == '#') continue;
			oceanMass.push_back(std::stod(tempBuffer));
		}
		f4.close();
	}
	catch (const std::exception& e)
	{
		f4.close();
		std::cout << "Initialisation of the 'oceanMass' vector failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	try
	{
		std::getline(f5, tempBuffer);
		if(tempBuffer[0] == '#') // first line may be a comment, starting with '#'
		{
			std::getline(f5, tempBuffer);
			model.setDt(stod(tempBuffer));
		}
		else
		{
			model.setDt(stod(tempBuffer));
		}
		std::getline(f5, tempBuffer);
		model.setThicknessSurfaceLayer(stod(tempBuffer));
		std::getline(f5, tempBuffer);
		model.setThicknessGasLayer(stod(tempBuffer));
		std::getline(f5, tempBuffer);
		model.setAtmMeanMolarMass(stod(tempBuffer));
		std::getline(f5, tempBuffer);
		model.setTemperatureGradient(stod(tempBuffer));
		std::getline(f5, tempBuffer);
		model.setSpecHeatCapacityAtm(stod(tempBuffer));
		std::getline(f5, tempBuffer);
		model.setPreviousCalcTime(stod(tempBuffer));

		std::getline(f5, tempBuffer);
		std::istringstream iss1(tempBuffer);
		std::string ss1;
		std::getline(iss1, ss1, ' ');
		model.setType(ss1);

		std::getline(f5, tempBuffer);
		std::istringstream iss2(tempBuffer);
		std::string ss2;
		std::getline(iss2, ss2, ' ');
		model.setFilename(ss2);

		f5.close();
	}
	catch (const std::exception& e)
	{
		f5.close();
		std::cout << "Initialisation of the model failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	try
	{
		std::getline(f6, tempBuffer);
		if(tempBuffer[0] == '#') // first line may be a comment, starting with '#'
		{
			std::getline(f6, tempBuffer);
			hoststar.setMass(stod(tempBuffer));
		}
		else
		{
			hoststar.setMass(stod(tempBuffer));
		}
		std::getline(f6, tempBuffer);
		hoststar.setRadius(stod(tempBuffer));
		std::getline(f6, tempBuffer);
		hoststar.setTemperature(stod(tempBuffer));

		f6.close();
	}
	catch (const std::exception& e)
	{
		f6.close();
		std::cout << "Initialisation of the hoststar failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	try
	{
		std::getline(f7, tempBuffer);
		if(tempBuffer[0] == '#') // first line may be a comment, starting with '#'
		{
			std::getline(f7, tempBuffer);
			exoplanet.setMass(stod(tempBuffer));
		}
		else
		{
			exoplanet.setMass(stod(tempBuffer));
		}
		std::getline(f7, tempBuffer);
		exoplanet.setRadius(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setAlbedo(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setEccentricity(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setSemiMajorAxis(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setDistance(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setOrbitalPeriod(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setPhaseAngle(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setNumberDensityAtSurface(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setVolumeRatioH2O(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setVolumeRatioCO2(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setVolumeRatioCH4(stod(tempBuffer));
		std::getline(f7, tempBuffer);
		exoplanet.setVolumeRatioN2O(stod(tempBuffer));

		f7.close();
	}
	catch (const std::exception& e)
	{
		f7.close();
		std::cout << "Initialisation of the exoplanet failed." << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;

}

int numberOfLines(std::ifstream& file)
{
	int i = 0;
	std::string line;
	while(getline(file, line)) i++;
	file.clear();
	file.seekg(0, std::ios::beg);
	return i;
}

void writeDataToFile(Model& model, Star& hoststar, Planet& exoplanet, std::string results, double diff)
{
	if(!model.getContinuation()) writeMetadata(model, hoststar, exoplanet, diff, ""); // Open the data file and write the (meta)data into it.
	std::fstream f(model.getFilename(), std::ios::app);
	if(!model.getContinuation()) f << "# iteration number | time [earth days] | orbital phase angle [rad] | surface air temperature [K] | surface temperature [K] | x position [m] | y position [m]" << std::endl;
	f << results;
	f.close();
}

int initVectors(Model& model,
		Star& hoststar,
		Planet& exoplanet,
		std::vector<double>& atmTemp,
		std::vector<double>& atmMass,
		std::vector<double>& oceanTemp,
		std::vector<double>& oceanMass,
		std::vector<double>& numDens,
		std::vector<double>& emitPower)
{
//	createTempVector(atmTemp, hoststar, exoplanet, model);
	atmTemp.reserve(model.getAtmTempSize());
	double h = hoststar.getTemperature() * std::sqrt(hoststar.getRadius() / exoplanet.getSemiMajorAxis()) * std::sqrt(std::sqrt((1 - exoplanet.getAlbedo()) / 4));
	atmTemp.push_back(h);
	for(int i = 1; i < model.getAtmTempSize(); i++)
	{
		h -= model.getTemperatureGradient() * model.getThicknessGasLayer();
		atmTemp.push_back(h);
	}

//	createNumDensVector(numDens, atmTemp, exoplanet, model);
	numDens.reserve(model.getAtmTempSize());
	for(int i = 0; i < (int)model.getAtmTempSize(); i++)
	{
		if(i == 0) numDens.push_back(exoplanet.getNumberDensityAtSurface());
		else numDens.push_back(atmNumberDensity(model, exoplanet, atmTemp[0], i));
	}

	// the number density is modified here for a better performance where the absorbed radiation is calculated
	for(unsigned int i = 0; i < numDens.size(); i++) numDens[i] *= ( exoplanet.getVolumeRatioH2O() * Constant::crossSectionH2O
                                                                   + exoplanet.getVolumeRatioCO2() * Constant::crossSectionCO2
                                                                   + exoplanet.getVolumeRatioCH4() * Constant::crossSectionCH4
                                                                   + exoplanet.getVolumeRatioN2O() * Constant::crossSectionN2O);

//	createMassVector(atmMass, atmTemp, numDens, exoplanet, model);
	atmMass.reserve(model.getAtmTempSize());
	atmMass.push_back(modelSurfaceMass(exoplanet, model));
	for(int i = 1; i < model.getAtmTempSize(); i++)
	{
		atmMass.push_back(massAtmLayer(model.getAtmMeanMolarMass(), numDens[i], i, exoplanet.getRadius(), model.getThicknessGasLayer()));
	}

	oceanTemp.reserve(model.getOceanTempSize());
	for(int i = 0; i < model.getOceanTempSize(); i++) oceanTemp.push_back(model.getOceanStartingTemp());

	oceanMass.reserve(model.getOceanTempSize());
	for(int i = 0; i < model.getOceanTempSize(); i++) oceanMass.push_back(modelSurfaceMass(exoplanet, model, i));

	emitPower.reserve(model.getAtmTempSize());
	emitPower.resize(model.getAtmTempSize());

	return 0;
}

void starPower(Model& model, Star& hoststar, Planet& exoplanet, double& temperature, double& mass)
{
	double specHeatCapacitySurface = getSpecHeatCapacitySurface(model.getType());

	double starPower = model.getDt() * (1 - exoplanet.getAlbedo()) * std::pow(exoplanet.getRadius()/exoplanet.getDistance(), 2)
			* Constant::stefanBoltzmann * M_PI * std::pow(hoststar.getRadius(), 2) * std::pow(hoststar.getTemperature(), 4)
			* (1 - hoststar.getIRPortion());

	changeTemperature(starPower, specHeatCapacitySurface, mass, temperature);
}
