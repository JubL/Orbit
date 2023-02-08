#ifndef MODELS_H
#define MODELS_H

#include "functions.h"
#include "obj.h"

#include <vector>
#include <mutex>

/**
	This method determines the black body temperatur for a body.
	The LaTeX representation of this formula is "T_\text{Oberfläche} &= T_\star \sqrt{\frac{R_\star}{d_\bullet}} \sqrt[4]{\frac{1-A}{4)}".
	@author Jubin Lirawi
	@param r_star The radius of the hoststar is given in meters.
	@param distance This is the distance between the hoststar and the (exo)planet; in meters.
	@param t_star The hoststars temperature is in kelvin.
	@param albedo The albedo of a body in a dimensionless quantity.
	@param epsilon Epsilon is the emissivity of a surface of a material. It is a dimensionless quantity. For black bodys it is equal to 1.
	@return The black body temperatur of a body is returned in kelvin.
*/
double temperature(double r_star, double distance, double t_star, double albedo, double epsilon = 1);

/**
	This function calculates the temperatur of an planet surface assuming the atmosphere is only one thin layer which reflects the infrared radiation back to the surface, but does not interact with the visible sunlight.
	The LaTeX representation of this formula is "T_\text{Surface} &= T_\star \sqrt{\frac{R_\star}{d_\bullet}} \sqrt[4]{\frac{1-A}{4(1 - \frac{\varepsilon}{2})}}".
	@author Jubin Lirawi
	@param r_star This is the radius of the star in meter.
	@param distance This is the distance between the star and the planet in meter.
	@param t_star The temperatur of the star is given in kelvin.
	@param albedo The albedo is the reflectivity of an surface and has no unit.
	@param epsilon Epsilon is the emissivity of a surface. It is dimensionless and has a value between 0 and 1. 1 is for black bodys (all radiation is absorbed) and 0 is for white bodys (all radiation is reflected).
	@return The surface temperatur is returned in kelvin.

*/
double oneLayerTemperature(double r_star, double distance, double t_star, double albedo, double epsilon);

/**
	This function is similar to the "oneLayerTemperatur" method, but instead of one atmospheric layer, it works with three of them.
	@author Jubin Lirawi
	@param r_star This is the radius of the star in meter.
	@param distance This is the distance between the star and the planet in meter.
	@param t_star The temperatur of the star is given in kelvin.
	@param albedo The albedo is the reflectivity of an surface and has no unit.
	@param epsilon Epsilon is the emissivity of a surface. It is dimensionless and has a value between 0 and 1. 1 is for black bodys (all radiation is absorbed) and 0 is for white bodys (all radiation is reflected).
	@return The surface temperatur is returned in kelvin.
*/
double threeLayerTemperature(double r_star, double distance, double t_star, double albedo, double epsilon);

/**
	This model calculates the surface temperature of a planet by considering the infrared emission of and the reflection of surface radiation by clouds.The presence of an atmosphere is neglected.
	The LaTeX representation of this formula is "T_\text{Surface} &= T_\star \sqrt{\frac{R_\star}{d_\bullet}} \sqrt[4]{\frac{1-A}{4\varepsilon(1-\nicefrac{c}{2})}}".
	@author Jubin Lirawi
	@param r_star This is the radius of the star in meter.
	@param distance This is the distance between the star and the planet in meter.
	@param t_star The temperatur of the star is given in kelvin.
	@param albedo The albedo is the reflectivity of an surface and has no unit.
	@param epsilon Epsilon is the emissivity of a surface. It is dimensionless and has a value between 0 and 1. 1 is for black bodys (all radiation is absorbed) and 0 is for white bodys (all radiation is reflected).
	@return The surface temperatur is returned in kelvin.
*/
double cloudyTemperature(double r_star, double distance, double t_star, double albedo, double epsilon, double cloudRatio = 0.6);

/**
	This model is to calculate the surface temperature of a planet.
	The atmosphere is devided in many thin layers, every which absorbs infrared radiation
	from the surface and heats up. The heat results in emission of radiation from which
	half goes down towards the surface. Clouds and winds are neglected.
	@author Jubin Lirawi
	@deprecated This method is deprecated.
	@param temperatureArray This vector holds the temperatures of all layers.
	@param layerMass This vector holds the masses of all layers.
	@param numDens This vector holds the number density of each layer.
	@param emitPower This vector holds the energy emited by an atmospheric layer.
	@param absPower This vector holds the energy absorbed by an atmospheric layer.
	@param model This object holds all the necessary model parameters.
	@param hoststar This is the Hoststar object.
	@param exoplanet This is the Exoplanet object.
	@param it_beg This is where the iteration shall start.
	@param it_end This is where the iteration shall end.
*/
[[deprecated("use 'greenhouse()' instead")]]
void manyLayers(std::vector<double>& temperatureArray,
		std::vector<double>& layerMass,
		std::vector<double>& numDens,
		Model& model,
		Star& hoststar,
		Planet& exoplanet,
		unsigned int it_beg = 0,
		unsigned int it_end = 0);

/*
	This model calculates the heat transport by conductivity within an ocean.
	@author Jubin Lirawi
	@param oceanTemp This vector holds the temperatures of all layers.
	@param oceanMass This vector holds the masses of all layers.
	@param temp This is the temperature of the surface layer. The surface temperature is determined by e.g. the 'manyLayers' model and another vector holds this value. Therefore 'temp0' has to be given to this method to match the surface temperature of that vector with that of this vector.
	@param exoplanet This is the planet object given to the method. From that the planetary radius is read out.
	@param model This object holds all the necessary model parameters.
*/
void deepOcean(std::vector<double>& oceanTemp, std::vector<double>& oceanMass, std::vector<double>& temp, Planet& exoplanet, Model& model);

/*
	This model calculates the heat transport by conductivity within the soil.
	@author Jubin Lirawi
	@param soilTemp This vector holds the temperatures of all layers.
	@param soilMass This vector holds the masses of all layers.
	@param exoplanet This is the planet object given to the method. From that the planetary radius is read out.
	@param model This object holds all the necessary model parameters.
*/
void soil(std::vector<double>& oceanTemp, std::vector<double>& oceanMass, Star& hoststar, Planet& exoplanet, Model& model);

/**
	This model is to account for the effect if the hoststar radiates in the infrared spectrum. In that case a portion of the radiation has to interact with the atmosphere before it reaches the surface.
	@author Jubin Lirawi
	@see manyLayers
	@param atmtemp This vector holds the temperatures of all layers.
	@param atmMass This vector holds the masses of all layers.
	@param numDens This vector holds the nmumber density of each layer.
	@param model This object holds all the necessary model parameters.
	@param hoststar This is the Hoststar object.
	@param exoplanet This is the Exoplanet object.
*/
void twoStreamApproxMod(std::vector<double>& atmTemp, std::vector<double>& atmMass, std::vector<double>& numDens, Model& model, Star& hoststar, Planet& exoplanet);

/**
	This model is an advancement of the 'manyLayers' model described above.
	@author Jubin Lirawi
	@param atmTemp This vector holds the temperatures of all layers.
	@param atmMass This vector holds the masses of all layers.
	@param numDens This vector holds the number density of each layer.
	@param emitPower This vector holds the energy emited by an atmospheric layer.
	@param model This object holds all the necessary model parameters.
	@param hoststar This is the Hoststar object.
	@param exoplanet This is the Exoplanet object.
	@param it_beg This is where the iteration shall start.
	@param it_end This is where the iteration shall end.
*/
void oldGreenhouse(std::vector<double>& atmTemp,
                   std::vector<double>& atmMass,
                   std::vector<double>& numDens,
                   std::vector<double>& emitPower,
                   Model& model,
                   Star& hoststar,
                   Planet& exoplanet,
                   unsigned int it_beg = 0,
                   unsigned int it_end = 0);

/**
	This model is an advancement of the 'manyLayers' model described above.
	@author Jubin Lirawi
	@param atmTemp This vector holds the temperatures of all atmospheric layers.
	@param atmMass This vector holds the masses of all atmopheric layers.
	@param soilTemp This vector holds the temperatures of all soil layers.
	@param soilMass This vector holds the masses of all soil layers.
	@param numDens This vector holds the number density of each layer.
	@param emitPower This vector holds the energy emited by an atmospheric layer.
	@param model This object holds all the necessary model parameters.
	@param exoplanet This is the Exoplanet object.
*/
void greenhouse(std::vector<double>& atmTemp,
                std::vector<double>& atmMass,
                std::vector<double>& soilTemp,
                std::vector<double>& soilMass,
                std::vector<double>& numDens,
                std::vector<double>& emitPower,
                Model& model,
                Planet& exoplanet);

#endif
