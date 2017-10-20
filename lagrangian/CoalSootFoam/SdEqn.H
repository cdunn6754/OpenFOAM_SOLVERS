{
    combustion->correct();

    fvScalarMatrix SdEqn
      (
       fvm::ddt(rho, Sd)
       + mvConvection->fvmDiv(phi, Sd)
       - fvm::laplacian(turbulence->muEff(), Sd)
       ==
       //       + combustion->R(Yi)
       + fvOptions(rho, Sd)

       );

    SdEqn.relax();

    fvOptions.constrain(SdEqn);

    //    YiEqn.solve(mesh.solver("Yi"));
    SdEqn.solve();

    fvOptions.correct(Yi);

    Sd.max(0.0);
}
