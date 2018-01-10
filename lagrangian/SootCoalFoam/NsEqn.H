{
  // volScalarField Sn
  //     (
  //         IOobject
  //         (
  //             "Sn",
  //             runTime.timeName(),
  //             mesh,
  //             IOobject::NO_READ,
  //             IOobject::NO_WRITE,
  //             false
  //         ),
  //         mesh,
  //         dimensionedScalar("zero", dimless/dimTime/dimVolume, 0.0)
  //     );


  
  fvScalarMatrix NsEqn
    (
     fvm::ddt(rho, Ns)
     + fvm::div(phi, Ns)
     - fvm::laplacian(turbulence->muEff()/SSc, Ns)
     ==
     sootSystem.sourceN(Ns)
     + fvOptions(rho, Ns)
     );

  NsEqn.relax();

  fvOptions.constrain(NsEqn);

  NsEqn.solve(mesh.solver("Ns"));

  fvOptions.correct(Ns);

  // Confusing but 
  // ensures minimum value of 0.0
  Ns.max(0.0); 

}
