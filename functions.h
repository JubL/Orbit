#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "models.h"
#include "constant.h"
#include "obj.h"
#include "pbar.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <mutex>

/**
	This function calculates an arithmetic mean value by summing all current values and dividing by the number of values.
	@param[out] meanValue The meanValue is given as an reference, so it can be manipulated during the execution of this function.
	@param currentValue This is one of many values to be averaged.
	@param n n is the number of values to be averaged and therefore a dimensionless quantity.
*/
void calcMean(double& meanValue, double currentValue, int n);

/**
	This function calculates the distance between a star and a planet on it's orbit
	The LaTeX representation for this formula is "r = \frac{a (1 - \varepsilon^2)}{1 + \varepsilon\cos\varphi}".
	@param e This is the numerical eccentricity of the orbit. It is dimensionless quantity and has therefore no unit.
	@param a This is the semi major axis of that orbit. It's given in meters.
	@param phi This is the phase angle between the planet on it's orbit and an arbitrary reference. Phi is given in radians.
	@return The distance is returned in meters.
*/
double distance(double e, double a, double phi);

/**
	This method calculates how much the orbital phase has changed after a certain amount of time has passed.
	The LaTeX representation for this formula is "\mathrm d\varphi = \frac{2\pi}{P}\;\frac{a^2 \; \sqrt{1 - \varepsilon^2}}{r^2} \; \mathrm dt".
	@param e This is the numerical eccentricity of the orbit. It is dimensionless quantity and has therefore no unit.
	@param a This is the semi major axis of that orbit. It's given in meters.
	@param r This is the distance between the hoststar and the exoplanet. It's given in meters.
	@param p This is the orbital period of the exoplanets orbit. It is given in the same unit as the time dt.
	@param dt This is the elapsed time for which the calculation is done. It is given in the same unit as the orbital period p.
	@return Returns the elapsed angle in respect to the elapsed time.
*/
double nextAngle(double e, double a, double r, double p, double dt);

/**
	This function calculates Wien's displacement law.
	@param temp This is the mean surface temperature of an object. It's unit is kelvin.
	@return The wavelenght of maximum intensity is returned in nanometers.
*/
double wiensDisplacementLaw(double temp);

/**
	This function calculates the radiative power of an object with a given temperature and surface area.
	@param area This is the surface area of an object witch is radiating. It is given in square meters.
	@param temperature The temperature of an object is put in in kelvin.
	@param epsilon Epsilon indicates whether an object behaves as an black, gray or white body. If epsilon is not given, it will be set to 1 forcing a black body behaviour. Epsilon is a dimensionless quantity.
	@return The return value is the power which is radiated by that object due to it's temperature. The unit is watt.
*/
double radiation(double area, double temperature, double epsilon = 1);

/**
	Using the 'radiation' method will lead to huge errors, especially if one multiplies that result with the time. This method calculates the radiation in small time steps and adjustes the temperature after each step.
	Not only is the radiative energy returned, also the temperature of that layer is calculated and set.
	@deprecated Do not use this function, as it is deprecated due to bad numerics. Use the 'radiateAndChangeTempRK4' function instead.
	@param[out] temperature The temperature of an object is put in in kelvin.
	@param area This is the surface area of an object witch is radiating. It is given in square meters.
	@param dt This is the time which shall elapse in total. It is given in seconds.
	@param layerMass This is the mass of a layer. It is given in kilogramms.
	@param specHeatCapacity The specific heat capacity is a material dependent constant. It's given in joule per kilogram and kelvin.
	@return The return value is the power which is radiated by that object due to it's temperature. The unit is watt.
*/
[[deprecated("use the RK4 equivalent instead")]]
double radiateAndChangeTemp(double& temperature, double area, double dt, double layerMass, double specHeatCapacity);

/**
	The 'radiateAndChangeTemp' method is an improvement of the 'radiation' method. This method improves the former by using a better ode solver, namely the RK4.
	The new temperature of an layer after radiating is calculated and set and moreover the radiated energy is returned.
	@param[out] temperature The temperature of an object is put in in kelvin.
	@param area This is the surface area of an object witch is radiating. It is given in square meters.
	@param dt This is the time which shall elapse in total. It is given in seconds.
	@param layerMass This is the mass of a layer. It is given in kilogramms.
	@param specHeatCapacity The specific heat capacity is a material dependent constant. It's given in joule per kilogram and kelvin.
	@return The return value is the power which is radiated by that object due to it's temperature. The unit is watt.
*/
double radiateAndChangeTempRK4(double& temperature, double area, double dt, double layerMass, double specHeatCapacity);

/**
	This function evaluates the new temperature of matter by adding the gained heat to the temperature it already has.
	@deprecated Do not use this function, as it is deprecated.
	@param heat This is the heat the matter has absorbed. It has the dimension joule.
	@param specHeatCapacity The specific heat capacity is a material dependent constant. It's given in joule per kilogram and kelvin.
	@param mass Mass has the unit kilogram.
	@param oldTemperature This is the temperature the object had had before, in kelvin.
	@return Returns the new temperature the object has now, in kelvin.
*/
[[deprecated("use 'changeTemperature()' instead")]]
double newTemperature(double heat, double specHeatCapacity, double mass, double oldTemperature);

/**
	This function evaluates and sets the new temperature of matter by adding the gained heat to the temperature it already has.
	@param heat This is the heat the matter has absorbed. It has the dimension joule.
	@param specHeatCapacity The specific heat capacity is a material dependent constant. It's given in joule per kilogram and kelvin.
	@param mass Mass has the unit kilogram.
	@param temperature This is the temperature the object had had before, in kelvin.
*/
void changeTemperature(double heat, double specHeatCapacity, double mass, double& temperature);

/**
	This function calculates how much power is absorbed by an layer of matter through which the power travels.
	@deprecated Do not use this function, as it is deprecated.
	@param initialQuantity This is the quantity to which the layer of matter is exposed.
	@param thickness The thickness of a layer through which a beam has to travel is given in meters.
	@param numberDensity1 The number density of absorbing matter type one, is the number of specified particles per volume and therefore the unit is the reciproc of cubic meters.
	@param crossSection1 The cross section is the effective area which indicates the likeliness of an absorbing effect. This value is for particles type one. It's given in square meters.
	@param numberDensity2 This is the number density of particles type two.
	@param crossSection2 This is the cross section of particles type two.
	@param numberDensity3 This is the number density of matter type three.
	@param crossSection3 This is the cross section of matter type three.
	@param numberDensity4 This is the number density of particles type four.
	@param crossSection4 This is the cross section of particles type four.
	@return The absorbed quantity is returned in the same unit as the initial quantity.
*/
[[deprecated("use 'absorb()' instead")]]
double absorbtion(double initialQuantity,
		double thickness,
		double numberDensity1,
		double crossSection1,
		double numberDensity2 = 0,
		double crossSection2 = 0,
		double numberDensity3 = 0,
		double crossSection3 = 0,
		double numberDenstiy4 = 0,
		double crossSection4 = 0);

/**
	This function calculates how much power is absorbed by an layer of matter through which the power travels.
	@param initialQuantity This is the quantity to which the layer of matter is exposed.
	@param thickness The thickness of a layer through which a beam has to travel is given in meters.
	@param numDens The number density of absorbing matter type multiplied with the volume ratio times the cross section. The unit is the reciproc of meters.
	@return The absorbed quantity is returned in the same unit as the initial quantity.
*/
double absorb(double initialQuantity,
		double thickness,
		double numDens);
/**
	This method returns the number density of an atmospheric gas for a given height.
	The LaTeX representation of this formula is "n(r) = n(r_\bullet) \left(\frac{r \; T(r_\bullet)}{r_\bullet \; \big(T(r_\bullet) + ar_\bullet-ar\big)}\right)^{- \frac{G \; m_\bullet \; M \; a}{R \; (T(r_\bullet) + ar_\bullet)^2}} \exp\left( - \frac{G \; m_\bullet \; M}{R \; \big(T(r_\bullet) + ar_\bullet\big)} \frac{r - r_\bullet}{r \; r_\bullet} \right)".
	@param model The model and hence it's parameters are given to this function.
	@param exoplanet The exoplanet and hence it's parameters are given to this function.
	@param surfaceTemperature This is the gas temperature at the surface. It's a boundary condition. Its dimension is kelvin.
	@param layerNumber This is the number of the considered layer. It's similar to the height above surface.
	@return The number density of the atmospheric gas is returned as the reciproc of cubic meters.
*/
double atmNumberDensity(Model& model, Planet& exoplanet, double surfaceTemperature, int layerNumber);

/*
	@param surfaceAtmPressure This is the gas pressure near the surface in pascal.
	@param surfaceAtmTemperature This is the gas temperature near the surface in kelvin.
	@return Returns the number density of the atmosphere near the surface as the reciproc of cubic meters.
*/
double surfaceAtmNumberDensity(double surfaceAtmPressure, double surfaceAtmTemperature);

/**
	This method calculates the mean molar mass of an atmospheric gas composition containing up to five different gases.
	If no values are given to this method, the atmospheric gas composition is assumed as follows: 78.08% Nitrogen, 20.95% Oxigen, 0.935328% water vapor, 0.0345% carbondioxide and 0.000172% Methan.
	@param volumeRatio1 This is the volume ratio of gas number 1. It is given as an dimensionless number.
	@param molarMass1 This is the molar mass of gas number 1. It is given in kg per mol.
	@param volumeRatio2 This is the volume ratio of gas number 2. It is given as an dimensionless number.
	@param molarMass2 This is the molar mass of gas number 2. It is given in kg per mol.
	@param volumeRatio3 This is the volume ratio of gas number 3. It is given as an dimensionless number.
	@param molarMass3 This is the molar mass of gas number 3. It is given in kg per mol.
	@param volumeRatio4 This is the volume ratio of gas number 4. It is given as an dimensionless number.
	@param molarMass4 This is the molar mass of gas number 4. It is given in kg per mol.
	@param volumeRatio5 This is the volume ratio of gas number 5. It is given as an dimensionless number.
	@param molarMass5 This is the molar mass of gas number 5. It is given in kg per mol.
	@return Returns the mean molar mass of a gas composition in kg per mol.
*/
double meanMolarMass(double volumeRatio1 = .7808,
			double molarMass1 =.028013,
			double volumeRatio2 = .2095,
			double molarMass2 = .031999,
			double volumeRatio3 = .00935328,
			double molarMass3 = .018015,
			double volumeRatio4 = .000345,
			double molarMass4 = .044010,
			double volumeRatio5 = .00000172,
			double molarMass5 = .016043);

/**
	This function is to calculate the mass of an thin, i.e. 1 meter thick, atmospheric layer.
	The LaTeX representation of this formula is "m &= M \; n \; \frac{1}{N_\text{A}} \; 4 \pi \; (r_\bullet + h)^2 \; \Delta r".
	@param molarMass Here we have the molar mass of a gas composition in kg per mol.
	@param numberDensity The number density of the atmospheric gas is needed to calculate the mass. It is in the reciproc of cubic meters.
	@param height This is the height in which the layer is.
	@param r_planet This is the radius of the considered planet in meters.
	@return Returns the mass of an atmospheric gas layer in kilograms.
*/
double massAtmLayer(double molarMass, double numberDensity, double height, double r_planet, double thickness);

/**
	This method recieves a vector, reserves the memory space and fills it with values.
	@deprecated Do not use this function, as it is deprecated.
	@param[out] temp This is the recieved and to be filled vector.
	@param hoststar A star object is given. It's needed for the calculations.
	@param exoplanet An exoplanet object is given. It's needed for the calculations.
	@param n This is the vectors lenght to be.
*/
[[deprecated("this is included in 'initVectors()'")]]
void createTempVector(std::vector<double>& temp, Star& hoststar, Planet& exoplanet, Model& model);

/**
	This method recieves a vector, reserves the memory space and fills it with values.
	@deprecated Do not use this function, as it is deprecated.
	@param[out] mass This is the recieved and to be filled vector.
	@param temp This is another vector needed for the calculations.
	@param exoplanet An exoplanet object is given. It's needed for the calculations.
	@param model The model and thus its parameters are given.
*/
[[deprecated("this is included in 'initVectors()'")]]
void createMassVector(std::vector<double>& mass, std::vector<double>& temp, std::vector<double>& numDens, Planet& exoplanet, Model& model);

/**
	This method recieves a vector, reserves the memory space and fills it with values.
	@deprecated Do not use this function, as it is deprecated.
	@param[out] mass This is the recieved and to be filled vector.
	@param temp This is another vector needed for the calculations.
	@param exoplanet An exoplanet object is given. It's needed for the calculations.
	@param model The model and thus its parameters are given.
*/
[[deprecated("this is included in 'initVectors()'")]]
void createNumDensVector(std::vector<double>& numDens, std::vector<double>& temp, Planet& exoplanet, Model& model);

/**
	This method is to calculate the mass of the surface of a planet. The thickness is considered to be one meter.
	@param exoplanet This is the planet object given to the method. From that the planetary radius is read out.
	@param type This string variable holds a keyword, so the method knows which density it has to use. Possible keywords are "water", "swamp", "earthlike", "rock", "desert" and "ice".
	@return Returns the mass of an one meter thick layer in kilogramms.
*/
double modelSurfaceMass(Planet& exoplanet, Model& model, int layer = 0);

/*
	This method returns the specific heat capacity of an desired material.
	@param type This is the type of the material. It is given as a string.
	@return Returns the specific heat capacity in Joule per kilogramm per Kelvin.
*/
double getSpecHeatCapacitySurface(std::string type);

/*
	This function calculates the energy which is transmitted by themal conductivity, id est the energy change of layer i.
	The LaTeX representation of this formula is "\mathrm dQ = -\lambda \; A \; \mathrm dt \; \frac{T_i - T_{i+1}}{d}".
	@param temp This is the vector which holds all the temperatures. It's needed to determine the differential temperature.
	@param i This is the number of the layer.
	@param exoplanet This is the planet object given to the method. From that the planetary radius is read out.
	@param model The model metadata is taken from here.
	@return The energy change (positiv or negative) is returned, so that the temperature change can be calculated as a next step.
*/
double calcThermalConduction(std::vector<double>& temp, int i, Planet& exoplanet, Model& model);

/**
	This function writes the metadata into the datafile.
	@param model The model metadata is taken from here.
	@param hoststar The metadata for the hoststar is taken from the hoststar object.
	@param exoplanet The metadata from the exoplanet is taken from the exoplant object.
	@param diff This is the time how long the calculation took.
	@param message This is an optional message.
*/
void writeMetadata(Model& model, Star& hoststar, Planet& exoplanet, double diff, std::string message = "");

/**
	This function is to write down all the value for the flux the planet receives from it's hoststar.
	@param starPower This is the powerfrom the hoststar received at the planet and corrected by the albedo.
	@param albedo This is the albedo, the reflection coefficient.
	@param planetRadius This is the radius of the planet.
*/
void noteStarFlux(double& starPower, double& albedo, double& planetRadius);

/**
	This method is to write out the neccessary data in case the calculation shall be continued.
	@param temp This is the vector containing the atmospheric layers temperatures.
	@param layerMass This is the vector containing the atmospheric layers masses.
	@param oceanTemp This is the vector containing the ocean layers temperatures.
	@param oceanMass This is the vector containing the ocean layers masses.
	@param model This object contains the model parameters like the timestep.
	@param hoststar This is the object representing the hoststar.
	@param exoplanet This is the object representing the exoplanet.
*/
void writeContinuationData(std::vector<double>& temp,
			std::vector<double>& layerMass,
			std::vector<double>& oceanTemp,
			std::vector<double>& oceanMass,
			Model& model,
			Star& hoststar,
			Planet& exoplanet);

/**
	This method is to load the neccessary data for the calculation to be continued.
	@param[out] temp This is the vector containing the atmospheric layers temperatures.
	@param[out] layerMass This is the vector containing the atmospheric layers masses.
	@param[out] oceanTemp This is the vector containing the ocean layers temperatures.
	@param[out] oceanMass This is the vector containing the ocean layers masses.
	@param[out] model This object contains the model parameters like the timestep.
	@param[out] hoststar This is the object representing the hoststar.
	@param[out] exoplanet This is the object representing the exoplanet.
	@return Returns 0 if successfull, -1 otherwise.
*/
int loadContinuationData(std::vector<double>& temp,
			std::vector<double>& layerMass,
			std::vector<double>& oceanTemp,
			std::vector<double>& oceanMass,
			Model& model,
			Star& hoststar,
			Planet& exoplanet);

/**
	This function returns the number of lines a text file has.
	@param file This is the filename of the to be counted text file.
	@retrun Returns the number of lines as an integer.
*/
int numberOfLines(std::ifstream& file);

/**
	This method writes the results out into a datafile. Furthermore some metadata is written out too.
	@param model This object contains the model parameters like the timestep.
	@param hoststar This is the object representing the hoststar.
	@param exoplanet This is the object representing the exoplanet.
	@param results This is the string containing the results.
	@param diff This is the time how long the calculation took.
*/
void writeDataToFile(Model& model, Star& hoststar, Planet& exoplanet, std::string results, double diff);

/**
	This method fills the vectors needed for the calculations and is used just for the readability of the main.
	Inside this method the number density vector is modified to save some time in the function which performs the absorbtion of the radiation.
	@param model This object contains the model parameters like the timestep.
	@param hoststar This is the object representing the hoststar.
	@param exoplanet This is the object representing the exoplanet.
	@param[out] atmTemp This is the vector containing the atmospheric layers temperatures.
	@param[out] atmMass This is the vector containing the atmospheric layers masses.
	@param[out] oceanTemp This is the vector containing the ocean layers temperatures.
	@param[out] oceanMass This is the vector containing the ocean layers masses.
	@param[out] numDens This vector contains the number density of each atmospheric layer.
	@param[out] emitPower This vector contains the energy emited by an atmospheric layer.
	@param[out] absPower This vector contains the energy absorbed by an atmospheric layer.
*/
int initVectors(Model& model,
		Star& hoststar,
		Planet& exoplanet,
		std::vector<double>& atmTemp,
		std::vector<double>& atmMass,
		std::vector<double>& oceanTemp,
		std::vector<double>& oceanMass,
		std::vector<double>& numDens,
		std::vector<double>& emitPower);
/**
	
	@param model This object contains the model parameters like the timestep.
	@param hoststar This is the object representing the hoststar.
	@param exoplanet This is the object representing the exoplanet.
	@param[out] temperature This is the temperature of the surface to be heated by the sun.
	@param mass This is the surfaces mass.
*/
void starPower(Model& model, Star& hoststar, Planet& exoplanet, double& temperature, double& mass);

#endif
