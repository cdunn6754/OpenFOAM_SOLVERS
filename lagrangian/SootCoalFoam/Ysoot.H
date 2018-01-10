{

    // Soot mass fraction
    volScalarField& Ysoot = Y[i];

    // // Soot mass fraction source 
    // volScalarField Sm
    //     (
    //         IOobject
    //         (
    //             "Sm",
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE,
    //             false
    //         ),
    //         mesh,
    //         dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    //     );

    // Soot schimdt number (need to look up again)
    const scalar SSc(60.0); // from phd thesis

    fvScalarMatrix YsootEqn
        (
            fvm::ddt(rho, Ysoot)
            + mvConvection->fvmDiv(phi, Ysoot)
            - fvm::laplacian(turbulence->muEff()/SSc, Ysoot)
            ==
            sootSystem.sourceY(Ysoot)
            + coalParcels.SYi(i, Ysoot)
            + combustion->R(Ysoot)
            + fvOptions(rho, Ysoot)
        );
    
    #include "NsEqn.H"

    // Adjust other species accordingly
    // C2H2.primitiveFieldRef()= C2H2
    //     - (R1*MwC2H2) * (runTime.deltaT()*mesh.V())
    //     - (R2*MwC2H2) * (runTime.deltaT()*mesh.V());
    // H2.primitiveFieldRef()= H2
    //     + (R2*MwH2*runTime.deltaT()*mesh.V())
    //     + (R2*MwH2*runTime.deltaT()*mesh.V());
    // O2.primitiveFieldRef()= O2
    //     - (0.5*R3*MwO2) * (runTime.deltaT()*mesh.V());
    // CO.primitiveFieldRef()= CO
    //     + (R3*MwCO) * (runTime.deltaT()*mesh.V())
    //     + (R4*MwCO) * (runTime.deltaT()*mesh.V());
    // OH.primitiveFieldRef()= OH
    //     - (R4*MwOH) * (runTime.deltaT()*mesh.V());
    // H.primitiveFieldRef() = H
    //     + (R4*MwH) * (runTime.deltaT()*mesh.V());
   

    YsootEqn.relax();

    fvOptions.constrain(YsootEqn);

    YsootEqn.solve(mesh.solver("Yi"));

    fvOptions.correct(Ysoot);

    Ysoot.max(0.0);
    Yt += Ysoot;
}
