Info<< "\nConstructing coal cloud" << endl;
DaemCoalCloud coalParcels
(
    "coalCloud1",
    rho,
    U,
    g,
    slgThermo
);

Info<< "\nConstructing limestone cloud" << endl;
basicThermoCloud limestoneParcels
(
    "limestoneCloud1",
    rho,
    U,
    g,
    slgThermo
);
