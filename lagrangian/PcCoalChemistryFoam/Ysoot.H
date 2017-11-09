{
    Info << "We are solving for SOOT now" << endl;

    // Soot mass fraction
    volScalarField& Ysoot = Y[i];

    // Molecular weights 
    const scalar& MwC2H2 = composition.W(composition.species()["C2H2"]);
    const scalar& MwSOOT = composition.W(composition.species()["SOOT"]);
    const scalar& MwH = composition.W(composition.species()["H"]);
    const scalar& MwO2 = composition.W(composition.species()["O2"]);
    const scalar& MwOH = composition.W(composition.species()["OH"]);
    const scalar& MwCO = composition.W(composition.species()["CO"]);
    const scalar& MwH2 = composition.W(composition.species()["H2"]);

    // Soot mass fraction source 
    volScalarField Sm
        (
            IOobject
            (
                "Sm",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
        );

    // Soot mass fraction reaction rate 1
    volScalarField R1
        (
            IOobject
            (
                "R1",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimMoles/dimTime/dimVolume, 0.0)
        );

    // Soot mass fraction reaction rate 2
    volScalarField R2
        (
            IOobject
            (
                "R2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimMoles/dimTime/dimVolume, 0.0)
        );

    // Soot mass fraction reaction rate 3
    volScalarField R3
        (
            IOobject
            (
                "R3",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimMoles/dimTime/dimVolume, 0.0)
        );

    // Soot mass fraction reaction rate 4
    volScalarField R4
        (
            IOobject
            (
                "R4",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimMoles/dimTime/dimVolume, 0.0)
        );

    // Soot Particle Diameter
    volScalarField dp
        (
            IOobject
            (
                "R4",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimLength, 0.0)
        );

    const dimensionedScalar NsSMALL("NsSmall", dimless/dimMass, SMALL);
    //Soot density
    const scalar rhos=1000;

    dp.primitiveFieldRef() = cbrt((6/pi) * (Ysoot/(rhos*(Ns+NsSMALL))));


    // Arrhenius Coefficients for reactions (from Kronenburg 2006)
    const dimensionedScalar Ta1("Ta", dimTemperature, 21000);
    const dimensionedScalar A1("A", dimless, 0.63e+04);
    const dimensionedScalar Ta2("Ta", dimTemperature, 12100);
    const dimensionedScalar A2("A", dimless, 0.75e+03);
    const dimensionedScalar Ta3("Ta", dimTemperature, 19800);
    const dimensionedScalar A3("A", dimless, 7.15e+02);
    const dimensionedScalar Ta4("Ta", dimTemperature, 0);
    const dimensionedScalar A4("A", dimless, 0.36);


    // Acetylene Mass Fraction
    volScalarField& C2H2(composition.Y("C2H2"));

    // Oxygen Mass Fraction
    volScalarField& O2(composition.Y("O2"));

    // Carbon Monoxide Mass Fraction
    volScalarField& CO(composition.Y("CO"));

    // Hydroxide Mass Fraction
    volScalarField& OH(composition.Y("OH"));

    // Radical Hydrogen Mass Fraction
    volScalarField& H(composition.Y("H"));

    // Hydrogen Mass Fraction
    volScalarField& H2(composition.Y("H2"));

    // Soot surface area per unit volume
    scalarField As(H2.size(), 0.0);
    As = pi*sqr(dp)*rho*Ns;

    // Calculate soot reaction rates
    R1.primitiveFieldRef() = (A1*exp(Ta1/T)) * ((C2H2*rho/MwC2H2));
    R2.primitiveFieldRef() = (A2*exp(Ta2/T)) * As * ((C2H2*rho/MwC2H2));
    R3.primitiveFieldRef() = (A3*exp(Ta3/T)*pow(T,0.5)) * As * ((O2*rho/MwO2));
    R4.primitiveFieldRef() = (A4*exp(Ta4/T)*pow(T,0.5)) * As * ((OH*rho/MwOH));

    // Assign the source term for soot
    Sm.primitiveFieldRef() = 2*R1*MwSOOT + 2*R2*MwSOOT
        - R3*MwSOOT - R4*MwSOOT;


    fvScalarMatrix YsootEqn
        (
            fvm::ddt(rho, Ysoot)
            + mvConvection->fvmDiv(phi, Ysoot)
            - fvm::laplacian(turbulence->muEff(), Ysoot)
            ==
            Sm
            + coalParcels.SYi(i, Ysoot)
            + combustion->R(Ysoot)
            + fvOptions(rho, Ysoot)
        );
    Info << "Defined soot equation, going to NS" << endl;

    #include "NsEqn.H"

    // Adjust other species accordingly
    C2H2.primitiveFieldRef()= C2H2
        - (R1*MwC2H2) * (runTime.deltaT()*mesh.V())
        - (R2*MwC2H2) * (runTime.deltaT()*mesh.V());
    H2.primitiveFieldRef()= H2
        + (R2*MwH2*runTime.deltaT()*mesh.V())
        + (R2*MwH2*runTime.deltaT()*mesh.V());
    O2.primitiveFieldRef()= O2
        - (0.5*R3*MwO2) * (runTime.deltaT()*mesh.V());
    CO.primitiveFieldRef()= CO
        + (R3*MwCO) * (runTime.deltaT()*mesh.V())
        + (R4*MwCO) * (runTime.deltaT()*mesh.V());
    OH.primitiveFieldRef()= OH
        - (R4*MwOH) * (runTime.deltaT()*mesh.V());
    H.primitiveFieldRef() = H
        + (R4*MwH) * (runTime.deltaT()*mesh.V());
    
            
    Info<<"Got through the soot" << endl;

    YsootEqn.relax();

    fvOptions.constrain(YsootEqn);

    YsootEqn.solve(mesh.solver("Yi"));

    fvOptions.correct(Ysoot);

    Ysoot.max(0.0);
    Yt += Ysoot;
}
