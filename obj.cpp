#include "obj.h"

Star::Star(double m, double r, double t, double ir)
{
	setMass(m);
	setRadius(r);
	setTemperature(t);
	setIRPortion(ir);
	std::cout << "Placing a star." << std::endl;
}

Sun::Sun()
{
	setMass(1);
	setRadius(1);
	setTemperature(5772);
	std::cout << "...the sun." << std::endl;
}

void Star::setMass(double m)
{
	if(m < 200) mass = m;
	else mass = m / Constant::massSun;
}

void Star::setRadius(double r)
{
	radius = r;
}

void Star::setTemperature(double t)
{
	temperature = t;
}

void Star::setIRPortion(double ir)
{
	if(ir <= 1 && ir >= 0) irPortion = ir;
	else throw std::invalid_argument("The stars ir flux portion must be between 0 and 1.");
}

double Star::getMass(int i)
{
	if(i == 0) return mass * Constant::massSun;
	else return mass;
}

double Star::getRadius(int i)
{
	if(i == 0) return radius * Constant::radiusSun;
	else return radius;
}

double Star::getTemperature()
{
	return temperature;
}

double Star::getIRPortion()
{
	return irPortion;
}

Planet::Planet(double m,
               double r,
               double a,
               double e,
               double sma,
               //double d,
               double period,
               double numDens,
               double ratioH2O,
               double ratioCO2,
               double ratioCH4,
               double ratioN2O)
{
	setMass(m);
	setRadius(r);
	setAlbedo(a);
	setEccentricity(e);
	setSemiMajorAxis(sma);
	setDistance(0);
	setOrbitalPeriod(period);
	setPhaseAngle(0);

	setNumberDensityAtSurface(numDens);
	setVolumeRatioH2O(ratioH2O);
	setVolumeRatioCO2(ratioCO2);
	setVolumeRatioCH4(ratioCH4);
	setVolumeRatioN2O(ratioN2O);
	std::cout << "Creating a planet." << std::endl;
}

Mercury::Mercury()
{
	setMass(0.0553);
	setRadius(0.383);
	setAlbedo(0.068);
    setEccentricity(0.20563069);
	setSemiMajorAxis(0.38709893);
	setDistance(0);
	setOrbitalPeriod(87.969);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Mercury." << std::endl;
}

Venus::Venus()
{
	setMass(0.815);
	setRadius(0.950);
	setAlbedo(0.77);
	setEccentricity(0.0067);
	setSemiMajorAxis(0.72333199);
	setDistance(0);
	setOrbitalPeriod(224.701);

	setNumberDensityAtSurface(9.0414 * pow(10, 23));
	setVolumeRatioH2O(0.00002);
	setVolumeRatioCO2(0.965);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Venus." << std::endl;
}

Earth::Earth()
{
	setMass(1);
	setRadius(1);
	setAlbedo(0.306);
	setEccentricity(0.01671022);
	setSemiMajorAxis(1.00000011);
	setDistance(0);
	setOrbitalPeriod(365.256);

	setNumberDensityAtSurface(2.867 * pow(10, 25));
	setVolumeRatioH2O(0.0003);
	setVolumeRatioCO2(0.00000345);
	setVolumeRatioCH4(0.0000000172);
	setVolumeRatioN2O(0.0000000031);
	std::cout << "...the Earth." << std::endl;
}

Moon::Moon()
{
	setMass(0.0123);
	setRadius(0.2727);
	setAlbedo(0.11);
	setEccentricity(0.01671022); // eccentricity of earth orbit!
	setSemiMajorAxis(1.00000011); // semimajoraxis of earth orbit!
	setDistance(0);
	setOrbitalPeriod(365.256);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...the moon." << std::endl;
}

Mars::Mars()
{
	setMass(0.107);
	setRadius(0.532);
	setAlbedo(0.25);
	setEccentricity(0.09341233);
	setSemiMajorAxis(1.52366231);
	setDistance(0);
	setOrbitalPeriod(686.980);

	setNumberDensityAtSurface(2.19 * pow(10, 23));
	setVolumeRatioH2O(0.0003);
	setVolumeRatioCO2(0.953);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Mars." << std::endl;
}

Ceres::Ceres()
{
	setMass(0.00015);
	setRadius(0.000157231335788);
	setAlbedo(0.09);
	setEccentricity(0.075823);
	setSemiMajorAxis(2.7675);
	setDistance(0);
	setOrbitalPeriod(1681.63);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...dwarfplanet Ceres." << std::endl;
}

Jupiter::Jupiter()
{
	setMass(317.83);
	setRadius(10.973);
	setAlbedo(0.343);
	setEccentricity(0.0489);
	setSemiMajorAxis(5.20336301);
	setDistance(0);
	setOrbitalPeriod(4332.589);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Jupiter." << std::endl;
}

Saturn::Saturn()
{
	setMass(95.16);
	setRadius(9.14);
	setAlbedo(0.342);
	setEccentricity(0.0565);
	setSemiMajorAxis(9.53707032);
	setDistance(0);
	setOrbitalPeriod(10759.22);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Saturn." << std::endl;
}

Uranus::Uranus()
{
	setMass(14.54);
	setRadius(3.981);
	setAlbedo(0.300);
	setEccentricity(0.0457);
	setSemiMajorAxis(19.19126393);
	setDistance(0);
	setOrbitalPeriod(30685.4);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Uranus." << std::endl;
}

Neptune::Neptune()
{
	setMass(17.15);
	setRadius(3.865);
	setAlbedo(0.290);
	setEccentricity(0.0113);
	setSemiMajorAxis(30.06896348);
	setDistance(0);
	setOrbitalPeriod(60189.0);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...Neptune." << std::endl;
}

Pluto::Pluto()
{
	setMass(0.0022);
	setRadius(0.186);
	setAlbedo(0.5);
	setEccentricity(0.2488);
	setSemiMajorAxis(39.48168677);
	setDistance(0);
	setOrbitalPeriod(90560.0);

	setVolumeRatioH2O(0);
	setVolumeRatioCO2(0);
	setVolumeRatioCH4(0);
	setVolumeRatioN2O(0);
	std::cout << "...dwarfplanet Pluto." << std::endl;
}

void Planet::setMass(double m)
{
	mass = m;
}

void Planet::setRadius(double r)
{
	radius = r;
}

void Planet::setAlbedo(double a)
{
	albedo = a;
}

void Planet::setEccentricity(double e)
{
	eccentricity = e;
}

void Planet::setSemiMajorAxis(double sma)
{
	semiMajorAxis = sma;
}

void Planet::setDistance()
{
	distance = Constant::astronomicalUnit * semiMajorAxis * (1 - pow(eccentricity, 2)) / (1 + eccentricity * std::cos(phaseAngle));
}

void Planet::setDistance(double d)
{
	distance = d;
}

void Planet::setOrbitalPeriod(double days)
{
	orbitalPeriod = days;
}

void Planet::setPhaseAngle(double phi)
{
	phaseAngle = phi;
}

void Planet::addPhaseAngle(double phi)
{
	phaseAngle += phi;
}

void Planet::calcPhaseAngle(double dt)
{
	phaseAngle += std::pow(Constant::astronomicalUnit * semiMajorAxis/distance, 2) * std::sqrt(1 - std::pow(eccentricity, 2)) * (2 * M_PI / (Constant::day * orbitalPeriod) ) * dt;
}

void Planet::setNumberDensityAtSurface(double value)
{
	numberDensityAtSurface = value;
}

void Planet::setVolumeRatioH2O(double value)
{
	volumeRatioH2O = value;
}

void Planet::setVolumeRatioCO2(double value)
{
	volumeRatioCO2 = value;
}

void Planet::setVolumeRatioCH4(double value)
{
	volumeRatioCH4 = value;
}

void Planet::setVolumeRatioN2O(double value)
{
	volumeRatioN2O = value;
}

double Planet::getMass(int i)
{
	if(i == 0) return mass * Constant::massEarth;
	else return mass;
}

double Planet::getRadius(int i)
{
	if(i == 0) return radius * Constant::radiusEarth;
	else return radius;
}

double Planet::getAlbedo()
{
	return albedo;
}

double Planet::getEccentricity()
{
	return eccentricity;
}

double Planet::getSemiMajorAxis(int i)
{
	if(i == 0) return semiMajorAxis * Constant::astronomicalUnit;
	else return semiMajorAxis;
}

double Planet::getDistance()
{
	if(distance == 0) return semiMajorAxis;
	else return distance;
}

double Planet::getOrbitalPeriod(int i)
{
	if(i == 0) return orbitalPeriod * Constant::day; // Could also be derived from Kepler III, but that means the Hoststar and the Exoplanet objects (their masses) have to be given to this method.
	else return orbitalPeriod;
}

double Planet::getPhaseAngle()
{
	return phaseAngle;
}

double Planet::getNumberDensityAtSurface()
{
	return numberDensityAtSurface;
}

double Planet::getVolumeRatioH2O()
{
	return volumeRatioH2O;
}

double Planet::getVolumeRatioCO2()
{
	return volumeRatioCO2;
}

double Planet::getVolumeRatioCH4()
{
	return volumeRatioCH4;
}

double Planet::getVolumeRatioN2O()
{
	return volumeRatioN2O;
}

Model::Model(unsigned int pDt,
             double pThicknessSurfaceLayer,
             double pThicknessGasLayer,
             double pAtmMeanMolarMass,
             std::string pType,
             double pTemperatureGradient,
             int pAtmTempSize,
             int pOceanTempSize,
             double pOceanStartingTemp,
             double pSpecHeatCapacityAtm,
             double pCalcTime,
             std::string pFilename,
             int pContinuation)
{
	dt = pDt;
	thicknessSurfaceLayer = pThicknessSurfaceLayer;
	thicknessGasLayer = pThicknessGasLayer;
	atmMeanMolarMass = pAtmMeanMolarMass;
	type = pType;
	temperatureGradient = pTemperatureGradient;
	atmTempSize= pAtmTempSize;
	oceanTempSize = pOceanTempSize;
	oceanStartingTemp = pOceanStartingTemp;
	specHeatCapacityAtm = pSpecHeatCapacityAtm;
	calcTime = pCalcTime;
	previousCalcTime = 0;
	filename = pFilename;
	continuation = pContinuation;

	std::cout << "Model parameters set." << std::endl;
}

void Model::setDt(unsigned int pDt)
{
	dt = pDt;
}

void Model::setThicknessSurfaceLayer(double pThicknessSurfaceLayer)
{
	thicknessSurfaceLayer = pThicknessSurfaceLayer;
}

void Model::setThicknessGasLayer(double pThicknessGasLayer)
{
	thicknessGasLayer = pThicknessGasLayer;
}

void Model::setAtmMeanMolarMass(double pAtmMeanMolarMass)
{
	atmMeanMolarMass = pAtmMeanMolarMass;
}

void Model::setType(std::string pType)
{
	type = pType;
}

void Model::setTemperatureGradient(double pTemperatureGradient)
{
	temperatureGradient = pTemperatureGradient;
}

void Model::setAtmTempSize(int pAtmTempSize)
{
	atmTempSize = pAtmTempSize;
}

void Model::setOceanTempSize(int pOceanTempSize)
{
	oceanTempSize = pOceanTempSize;
}

void Model::setOceanStartingTemp(double pOceanStartingTemp)
{
	oceanStartingTemp = pOceanStartingTemp;
}

void Model::setSpecHeatCapacityAtm(double pSpecHeatCapacityAtm)
{
	specHeatCapacityAtm = pSpecHeatCapacityAtm;
}

void Model::setCalcTime(double pCalcTime) // in local years
{
	calcTime = pCalcTime;
}

void Model::setPreviousCalcTime(double pPreviousCalcTime) // in earth years
{
	previousCalcTime = pPreviousCalcTime;
}

void Model::setFilename(std::string pFilename)
{
	filename = pFilename;
}

void Model::setContinuation(int pContinuation) // 0 = false; 1 = true
{
	continuation = pContinuation;
}

unsigned int Model::getDt()
{
	return dt;
}

double Model::getThicknessSurfaceLayer()
{
	return thicknessSurfaceLayer;
}

double Model::getThicknessGasLayer()
{
	return thicknessGasLayer;
}

double Model::getAtmMeanMolarMass()
{
	return atmMeanMolarMass;
}

std::string Model::getType()
{
	return type;
}

double Model::getTemperatureGradient()
{
	return temperatureGradient;
}

int Model::getAtmTempSize()
{
	return atmTempSize;
}

int Model::getOceanTempSize()
{
	return oceanTempSize;
}

double Model::getOceanStartingTemp()
{
	return oceanStartingTemp;
}

double Model::getSpecHeatCapacityAtm()
{
	return specHeatCapacityAtm;
}

double Model::getCalcTime() // in local years
{
	return calcTime;
}

double Model::getPreviousCalcTime() // in earth years
{
	return previousCalcTime;
}

std::string Model::getFilename()
{
	return filename;
}

int Model::getContinuation() // 0 = false; 1 = true
{
	return continuation;
}

