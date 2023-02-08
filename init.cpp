Sun hoststar;
Earth exoplanet;

exoplanet.setEccentricity(.0);

Model model;

model.setDt(20);
model.setCalcTime(1);
model.setAtmTempSize(1);
//model.setThicknessGasLayer(110)
//model.setThicknessSurfaceLayer(.05);
//model.setOceanTempSize(200);
model.setOceanStartingTemp(254.05); // <-- blackbody temperature of the earth
model.setFilename("test.dat");
