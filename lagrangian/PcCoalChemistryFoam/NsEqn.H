{
  volScalarField Sn
      (
          IOobject
          (
              "Sn",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false
          ),
          mesh,
          dimensionedScalar("zero", dimless/dimTime/dimVolume, 0.0)
      );

  // Avagadro's number
  const scalar Na(6.232e26);
  // Steffan Boltzmann Constant
  const scalar sigB(1.38e-23);
  // Minumum carbon atoms for soot particle nucleation
  const scalar nc(60);
  // Agglomeration constant (constants from Kronenburg 2006)
  const scalar Ca(9.0);


  // Source for number of soot particles
  Sn.primitiveFieldRef() = 2 * R1.primitiveField() * (Na/nc) 
      - 2*Ca*dp.primitiveField()*pow(6*sigB*T.primitiveField()/rhos,0.5) 
      * pow(rho*Ns.primitiveField(),2);
  
  fvScalarMatrix NsEqn
    (
     fvm::ddt(rho, Ns)
     + fvm::div(phi, Ns)
     - fvm::laplacian(turbulence->muEff()/SSc, Ns)
     ==
     Sn
     + fvOptions(rho, Ns)
     );

  NsEqn.relax();

  fvOptions.constrain(NsEqn);

  NsEqn.solve(mesh.solver("Ns"));

  fvOptions.correct(Ns);

}
