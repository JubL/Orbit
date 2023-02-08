#ifndef OBJ_H
#define OBJ_H

#include "constant.h"

#include <stdexcept>
#include <iostream>
#include <cmath>

/**
	This class represents stars in general.
*/
class Star
{
	public:
		Star(double m = 1,
             double r = 1,
             double t = 5772,
             double ir = 0);
		void setMass(double m); /** set in solar masses or kilogram */
		void setRadius(double r); /** set in solar radii */
		void setTemperature(double t); /** set in kelvin */
		void setIRPortion(double ir); /** set as dimensionless number between 0 and 1 */
		double getMass(int i = 0); /** in kg (default) or solar masses */
		double getRadius(int i = 0); /** in meter (default) or solar radii */
		double getTemperature(); /** in Kelvin */
		double getIRPortion(); /** dimensionless between 0 and 1 */
	protected:
		double mass; /** in solar masses */
		double radius; /** in solar radii */
		double temperature; /** in Kelvin */
		double irPortion; /** dimensionless between 0 and 1 */
};

/**
	This class represents a special star: our sun. It is derived from the star class.
*/
class Sun : public Star
{
	public:
		Sun();
};

/**
	This class represents planets in general.
*/
class Planet
{
	public:
		Planet(double m = 1,
               double r = 1,
               double a = 0.306,
               double e = 0.01671022,
               double sma = 1.00000011,
               //double d = 0,
               double period = 365.256,
               double numDens = 2.687 * pow(10, 25),
               double ratioH2O = 0.0003,
               double ratioCO2 = 0.00000345,
               double ratioCH4 = 0.0000000172,
               double ratioN2O = 0.0000000031);
		void setMass(double m); /** set in earth masses */
		void setRadius(double r); /** set in earth radii */
		void setAlbedo(double a); /** as dimensionless number */
                void setEccentricity(double e); /** as dimensionless number */
		void setSemiMajorAxis(double sma); /** set in astronomical units */
		void setDistance();
		void setDistance(double d); /** set in meters */
		void setOrbitalPeriod(double days); /**set in earth days */
		void setPhaseAngle(double phi); /** as dimensionless number */
		void addPhaseAngle(double phi); /** as dimensionless number */
		void calcPhaseAngle(double dt); /** in seconds */
		void setNumberDensityAtSurface(double value); /** in reciproc cubic meters */
		void setVolumeRatioH2O(double value); /** as dimensionless number */
		void setVolumeRatioCO2(double value); /** as dimensionless number */
		void setVolumeRatioCH4(double value); /** as dimensionless number */
		void setVolumeRatioN2O(double value); /** as dimensionless number */
		double getMass(int i = 0); /** in kg (default) or earth masses */
		double getRadius(int i = 0); /** in meter (default) or earth radii */
		double getAlbedo(); /** as dimensionless number */
		double getEccentricity();  /** as dimensionless number */
		double getSemiMajorAxis(int i = 0); /** in meters (default) or in astronomical units */
		double getDistance(); /** in meters */
		double getOrbitalPeriod(int i = 0); /** in seconds (default) or earth days */
		double getPhaseAngle(); /** as dimensionless number */
		double getNumberDensityAtSurface();
		double getVolumeRatioH2O(); /** as dimensionless number */
		double getVolumeRatioCO2(); /** as dimensionless number */
		double getVolumeRatioCH4(); /** as dimensionless number */
		double getVolumeRatioN2O(); /** as dimensionless number */
	protected:
		double mass; /** in earth masses */
		double radius; /** in earth radii */
		double albedo; /** as dimensionless number */
		double eccentricity; /** numerical eccentricity */
		double semiMajorAxis; /** in astronomical units */
		double distance; /** in meter */
		double orbitalPeriod; /** in earth days */
		double phaseAngle; /** as dimensionless number */
		double numberDensityAtSurface; /** in reciproc cubic meters */
		double volumeRatioH2O; /** as dimensionless number */
		double volumeRatioCO2; /** as dimensionless number */
		double volumeRatioCH4; /** as dimensionless number */
		double volumeRatioN2O; /** as dimensionless number */
};

/**
	This class represents a special planet: Mercury
*/
class Mercury : public Planet
{
	public:
		Mercury();
};

/**
	This class represents a special planet: Venus
*/
class Venus : public Planet
{
	public:
		Venus();
};

/**
	This class represents a special planet: Earth
*/
class Earth : public Planet
{
	public:
		Earth();
};

/**
	This class represents a special planet: our Moon
*/
class Moon : public Planet
{
	public:
		Moon();
};

/**
	This class represents a special planet: Mars
*/
class Mars : public Planet
{
	public:
		Mars();
};

/**
	This class represents a special celestial body: dwarf planet Ceres
*/
class Ceres : public Planet
{
	public:
		Ceres();
};

/**
	This class represents a special planet: Jupiter
*/
class Jupiter : public Planet
{
	public:
		Jupiter();
};

/**
	This class represents a special planet: Saturn
*/
class Saturn : public Planet
{
	public:
		Saturn();
};

/**
	This class represents a special planet: Uranus
*/
class Uranus : public Planet
{
	public:
		Uranus();
};

/**
	This class represents a special planet: Neptun
*/
class Neptune : public Planet
{
	public:
		Neptune();
};

/**
	This class represents a special celestial body: dwarf planet Pluto
*/
class Pluto : public Planet
{
	public:
		Pluto();
};

class Model
{
	public:
		Model(unsigned int pDt = 200,
              double pThicknessSurfaceLayer = .05,
              double pThicknessGasLayer = 110,
              double pAtmMeanMolarMass = .028586,
              std::string pType = "water",
              double pTemperatureGradient = .0065,
              int pAtmTempSize = 100,
              int pOceanTempSize = 200,
              double pOceanStartingTemp = 288,
              double pSpecHeatCapacityAtm = 1005,//717.6,
              double PCalcTime = 2.0,
              std::string pFilename = "test.dat",
              int pContinuation = 0);
		void setDt(unsigned int pDt);
		void setThicknessSurfaceLayer(double pThicknessSurfaceLayer);
		void setThicknessGasLayer(double pThicknessGasLayer);
		void setAtmMeanMolarMass(double pAtmMeanMolarMass);
		void setType(std::string pType);
		void setTemperatureGradient(double pTemperatureGradient);
		void setAtmTempSize(int pAtmTempSize);
		void setOceanTempSize(int pOceanTempSize);
		void setOceanStartingTemp(double pOceanStartingTemp);
		void setSpecHeatCapacityAtm(double pSpecHeatCapacityAtm);
		void setCalcTime(double pCalcTime); /** set in planetary years, i.e. revolutions around the hoststar */
		void setPreviousCalcTime(double pPreviousCalcTime); /* in earth years */
		void setFilename(std::string pFilename);
		void setContinuation(int pContinuation); /** 0 = false; 1 = true */
		unsigned int getDt();
		double getThicknessSurfaceLayer();
		double getThicknessGasLayer();
		double getAtmMeanMolarMass();
		std::string getType();
		double getTemperatureGradient();
		int getAtmTempSize();
		int getOceanTempSize();
		double getOceanStartingTemp();
		double getSpecHeatCapacityAtm();
		double getCalcTime(); /** in earth years */
		double getPreviousCalcTime(); /** in earth years */
		std::string getFilename();
		int getContinuation(); /** 0 = false; 1 = true */
	private:
		unsigned int dt;
		double thicknessSurfaceLayer;
		double thicknessGasLayer;
		double atmMeanMolarMass;
		std::string type;
		double temperatureGradient;
		int atmTempSize;
		int oceanTempSize;
		double oceanStartingTemp;
		double specHeatCapacityAtm;
		double calcTime;
		double previousCalcTime;
		std::string filename;
		int continuation;
};

#endif
