01-08-18:

PCCoalChemistryFoam:

PcCoalChemistryFoam uses the PCReactingMultiphase library to permit the use of PC Coal Lab devolatilization parameters and its general devolatilization scheme. In PcCoalChemistryFoam soot and tar are produced (along with some hydrocarbons and other species) by primary devolatilization of the coal parcels. The primary products and the soot are then released into the gas phase (we are modeling soot as just a gas species and not lagrangian in any way). The soot was previously just transported passively. 

In the latest version of PcCoalChemistryFoam the two equations from Kronenburg CMC paper for soot mass fraction and particle number density evolution are incorporated. The source terms from those equations depend on other species concentrations as well as other thermo properties like temperature. The two transport equations are solved as the other species are but the solution their source terms are segregated from the overall chemistry evaluation. That is to say that the all other chemical species concentrations are frozen for the evaluation of the soot equation sources instead of being coupled. In PcCoalChemistryFoam the soot source equations are solved explicitly (with a very basic forward euler method). This led to some problems with overshooting a reasonable mass fraction for soot (i.e. anything outside of [0,1] is not going to work well). 

With this explicit method for evaluating the source terms it is easy to determine the concomitant change in the other species invovled (e.g. C2H2 is used for soot nucleation so we find that rate while doing the explicit step anyways). We then add that 'other specie' source to the species field directly. In the future we will add it into the source side of the equation instead.

The tar produced in primary devolatilization remains with the coal parcel that produced it (obviously a potentially large modeling assumption to make). This let us use the PC coal lab secondary devolatilization/pyrolysis parameters to govern the breakdown of that resident tar into smaller species (hydrocarbons and CO, CO2, ...) that actually make sense to store in the gas phase. As the parcel resident tars breakdown (again these tars are modeled as staying 'within' the parcel) their breakdown products are released to the gas phase as in primary devolatilization. 

SootCoalFoam:

Built on PCCoalChemistryFoam. All of the devolatilization stuff is the same. Now we are trying to solve for the source terms of the soot equations in an implicit manner. The source terms are functions of eachother, the thermo state, some of the other species concentrations and constants. For example in the case of the soot mass fraction 'Y_soot' we have something like dY_soot/dt = f(Y_soot,N_soot,Temp, Pressure, density, Y_i, C) where N_soot is the soot particle number density, Y_i are various other species concentrations and C is a constant (or list of)). 

We solve these two coupled equations by freezing all other variables (species and thermo). They are still non-linear so we use a pretty basic newton's method to solve them (using the overall gas phase time step). Once the next time step values are determined we can find the source over the integration interval by subtracting the previous value and dividing by the integration time step. Then, for every cell (all of this is done in an individual cell by the way), we just add the source term into the two soot equations. 

This implicit method approach makes it difficult to determine an appropriate change/source in the other chemical species. Particularly since each individual soot source equation has multiple terms it is not easy to assign the change in the soot variable to a particular term. 

One option would be to just use the explicit method as in PCCoalChemistryFoam. This is undesirable since the whole point of using implicit methods is to avoid any possible overshoot on the other specie predictions. The other obvious extreme is to include these species as varibles in the implicit method and make it a 5x5 system or whatever rather than a 2x2 system with the others frozen. 
