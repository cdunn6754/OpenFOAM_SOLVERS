
{

  // added clint 11-17 from brown-fletcher 1998
  scalar coalSchmidt = 700;

  fvScalarMatrix NsEqn
    (
     fvm::ddt(rho, Ns)
     + fvm::div(phi, Ns)
     - fvm::laplacian(turbulence->muEff()/coalSchmidt, Ns)
     ==
     coalParcels.SNs(Ns)
     + combustion->R(Ns)
     + fvOptions(rho, Ns)
     );

  NsEqn.relax();

  fvOptions.constrain(NsEqn);

  NsEqn.solve(mesh.solver("Ns"));

  fvOptions.correct(Ns);

}
