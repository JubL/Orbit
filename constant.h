#ifndef CONSTANT_H
#define CONSTANT_H

#include <cmath>

/**
 * This struct provides fundamental physical and natural constants.
 */
struct Constant
{
    static constexpr double stefanBoltzmann       = 5.670367e-8;
    static constexpr double gasConstant           = 8.3144598;
    static constexpr double avogadroConstant      = 6.022140857e23;
    static constexpr double boltzmannConstant     = 1.38064852e-23;
    static constexpr double gravitationalConstant = 6.67408e-11;              /** This is the gravitational constant in cubic meters per kilogram per square second. https://physics.nist.gov/cgi-bin/cuu/Value?bg */
    static constexpr double astronomicalUnit      = 149.5978707e9;            /** This is the Astronomic Unit in meters. https://cneos.jpl.nasa.gov/glossary/au.html */
    static constexpr double radiusSun             = 695700000;                /** This is the radius of the Sun in meters. https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html */
    static constexpr double radiusJupiter         = 69911000;                 /** This is the radius of Jupiter in meters. https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html */
    static constexpr double radiusEarth           = 6371000;                  /** This is the volumetric mean radius of the Earth in meters. https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html */
    static constexpr double radiusMoon            = 1737400;                  /** This is the volumetric mean radius of the Moon in meters. https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html */
    static constexpr double massSun               = 1.9885e30;                /** This is the mass of the Sun in kilogramms. https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html */
    static constexpr double massJupiter           = 1898.19e24;               /** This is the mass of Jupiter in kilogramms. https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html */
    static constexpr double massEarth             = 5.9723e24;                /** This is the mass of the Earth in kilogramms. https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html */
    static constexpr double massMoon              = 0.07346e24;               /** This is the mass of the Moon in kilogramms. https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html */
    static constexpr double day                   = 0620 * 0330;              /** This is the length of one earth day in seconds. */
    static constexpr double year                  = 0x16D.4189374BCP0L * day; /** This is the length of one earth year in seconds. */
    static constexpr double crossSectionH2O       = 2.22547e-25;              /** from exomol.com */
    static constexpr double crossSectionCO2       = 1.14906e-24;
    static constexpr double crossSectionCH4       = 1.83999e-25;
    static constexpr double crossSectionN2O       = 0;
private:
    Constant() = delete;
};

#endif
