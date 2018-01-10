/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

//using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <math.h>
class EulerImplicitSystem
{

    // true constants
    const scalar rho_s = 2000.; // soot density [kg/m^3]
    const scalar Na = 6.232e26; //Avagadros number [particles/kMol]
    const scalar sigma = 1.38e-23; // boltzman [J/K]
    const scalar Ca = 9.; //agglomeration constant [-]
    const scalar n_c = 60.; // minimum particle number for agglomeration [particles]?
    const scalar pi = 3.1415;

    // make the proper dimension sets for the sources
    dimensionSet Y_source_dims; //[kg/(m^3 * s)]
    dimensionSet N_source_dims; // [1/(m^3 * s)]

    // molecular weights
    const scalar MW_soot; // Molecular weight of soot [kg/kmol]
    const scalar MW_c2h2;
    const scalar MW_oh;
    const scalar MW_o2;

    // Thermo fields
    const scalarField& rho;
    const scalarField& T;
    const scalarField& Y_c2h2;
    const scalarField& Y_oh;
    const scalarField& Y_o2;
    
    // Soot fields
    const scalarField& Y_s;
    const scalarField& N_s;

    // Soot source fields for main Yeqn
    scalarField Y_source;
    scalarField N_source;
    

    // Constants are determined from examination of 
    // soot equation source equations while
    // only considering Y_s and N_s to be varible

    // Particle surface area constant
    scalar C_A;

    // Constants for the equations form of 'C_equation_term'
    scalar C_1_1;
    scalar C_1_2;
    scalar C_1_3;
    scalar C_1_4;
    
    scalar C_2_1;
    scalar C_2_2;

    // Set the constants based on new time step species/rho/T
    void setConstants(const label& cellNumber) 
        {
            // Get quantities for particular cell
            const scalar& rho_ = rho[cellNumber];
            const scalar& T_ = T[cellNumber];
            const scalar& Y_c2h2_ = Y_c2h2[cellNumber];
            const scalar& Y_oh_ = Y_oh[cellNumber];
            const scalar& Y_o2_ = Y_o2[cellNumber];
            
            // Particle surface area constant, used in both equations
            C_A = rho_ * pi * Foam::pow(6/(pi * rho_s), (2./3.));

            // Soot mass fraction equation (Y*rho_/MW) = concentrations [kmol/m^3]
            C_1_1 = 2. * (6.3e3* Foam::exp(-21000./T_)) * 
                (Y_c2h2_*rho_/MW_c2h2) * MW_soot;
            C_1_2 = 2. * (7.5e2 * Foam::exp(-12100./T_)) * C_A * 
                (Y_c2h2_*rho_/MW_c2h2) * MW_soot;
            C_1_3 = (7.15e2 * Foam::sqrt(T_) * Foam::exp(-19800./T_)) * 
                C_A * (Y_o2_*rho_/MW_o2) * MW_soot;
            C_1_4 = 0.36 * Foam::sqrt(T_) * C_A * 
                (Y_oh_*rho_/MW_oh) * MW_soot;

            // Soot particle number density equation
            C_2_1 = 2. * (6.3e3 * Foam::exp(-21000./T_)) *
                (Y_c2h2_*rho_/MW_c2h2) * (Na/n_c);
            C_2_2 = 2. * Ca * Foam::pow(6./(pi * rho_s), (1./6.))
                * Foam::sqrt((6. * sigma * T_)/rho_s) 
                * Foam::sqr(rho_);
        }


    // Member Functions
        
    label nEqns() const
        {
            return 2;
        }

    // The value of the derivatives given Y and N
    // I.E. dY/dt = f(Y,N) with f taken from Kronenburg 
    // CMC paper
    void derivatives
    (
        const scalar& Y, 
        const scalar& N,
        scalarField& derivative
    )
        {
            
            // Evolution of soot mass fraction y[0]
            derivative[0] = (C_1_1 + 
            (C_1_2 - C_1_3 - C_1_4) * (Foam::pow(Y,(2./3.)) * Foam::pow(N,(1./3.))));

            // Evolution of soot number density (particles/m^3)
            derivative[1] = (C_2_1 - 
            C_2_2 * (Foam::pow(Y,(1./6.)) * Foam::pow(N,(11./6.))));
        }

    // Evaluate the newton equations for the euler implicit form of the 
    // two source equations i.e. 'Y^i - Y^(i-1) - dt(dY/dt^i) = 0
    // for a single cell.
    // Assumes
    // Y_i and N_i hold the values at the integration start time.
    // 
    // Y,N - values for Y_soot and N_soot (guessed values most often)
    // dt - Euler implicit integration time step
    // fY, fN - values of the newton method equations evaluated at Y and N
    void newtonEquations
    (
        const scalar& Y,
        const scalar& N,
        const scalar& Y_i,
        const scalar& N_i,
        const scalar dt,
        scalar& fY,
        scalar& fN
    ) const
        {
            // Evolution of soot mass fraction y[0] equation
            fY = Y - Y_i - dt * (C_1_1 + 
            (C_1_2 - C_1_3 - C_1_4) * (Foam::pow(Y,(2./3.)) * Foam::pow(N,(1./3.))));

            // Evolution of soot number density (particles/m^3) equation
            fN = N - N_i - dt * (C_2_1 - 
            C_2_2 * (Foam::pow(Y,(1./6.)) * Foam::pow(N,(11./6.))));

        }

    // Find the Jacobian of the newton method system
    // i.e. the derivatives of the functions from
    // newtonEquations
    // again this peforms the calculation for a single cell.
    void jacobian
    (
        const scalar& Y,
        const scalar& N,
        const scalar dt,
        scalarSquareMatrix& J,
        scalarSquareMatrix& invJ
    ) const
        {   

        // Versions to be used in denominators
        const scalar Ysafe = Foam::max(Y + Foam::SMALL, Foam::SMALL);        
        const scalar Nsafe = Foam::max(N + Foam::SMALL, Foam::SMALL);

        // set the jacobian entries
        J(0, 0) = 1 - dt*(C_1_2 - C_1_3 - C_1_4) *
            (2./3.) * Foam::pow(Ysafe,(-1./3.)) * Foam::pow(N,(1./3.));

        J(0, 1) = - dt*(C_1_2 - C_1_3 - C_1_4) *
            Foam::pow(Y,(2./3.)) * (1./3.) * Foam::pow(Nsafe,(-2./3.));

        J(1, 0) = - dt * C_2_2 * 
           (1./6.) * Foam::pow(Ysafe,(-5./6.)) * Foam::pow(N,(11./6.));

        J(1, 1) = 1 - dt * C_2_2 * 
            Foam::pow(Y,(1./6.)) * (11./6.) * Foam::pow(Nsafe,(5./6.));

         
        // now invert the jacobian
        this->invertJacobian(J,invJ);

    }

    void invertJacobian
    (
        const scalarSquareMatrix& J,
        scalarSquareMatrix& invJ
    ) const
        {
            // get inverse of 2x2 jacobian with adjugate/determinant method
            scalar detJ = J(0,0)*J(1,1) - J(0,1)*J(1,0);

            if (fabs(detJ) < Foam::SMALL)
            {
                Info << "Small determinant: " << detJ << endl;
                detJ = Foam::SMALL;
            }

            invJ(0,0) = J(1,1)/detJ;
            invJ(0,1) = -J(0,1)/detJ;
            invJ(1,0) = -J(1,0)/detJ;
            invJ(1,1) = J(0,0)/detJ;
            
        }

    void explicitEuler
    (
        scalar& Y_i,
        scalar& N_i,
        const scalar& dt
    )
        {
            // get the derivatives at this location
            scalarField derivative(this->nEqns());
            this->derivatives(Y_i,N_i,derivative);

            // Take explicit euler step
            Y_i = Y_i + derivative[0]*dt;
            N_i = N_i + derivative[1]*dt;


            // Force them to be positive
            if (Y_i <= 0. || N_i <= 0.)
            {
                Y_i = 0.0;
                N_i = 0.0;
            }
                
        }

    // Integrate over timestep dt.
    // Once Y^(t+dt),N ^(t + dt) are found use them to find the 
    // source terms. Update sources.
    void correct
    (
        const scalar dt,
        const scalar relTol = 1e-6
    )
        {

            // Counter for number of cells in which
            // the newton method fails to converge
            label newton_method_failures(0);

            // loop through all cells
            forAll(Y_s, cellNumber)
            {
                //Info << "\n\nNew cell" << endl;

                // Reset the values of 'constants'
                // for this cell
                this->setConstants(cellNumber);

                // set initial values for this cell (to remain const)
                const scalar Y_0 = Y_s[cellNumber];
                const scalar N_0 = N_s[cellNumber];
                
                // set intial values for this cell
                scalar Y_i = Y_s[cellNumber];
                scalar N_i = N_s[cellNumber];

                // Info << "Y_next, N_next: " << Y_i
                //     << ", " << N_i << endl;

                // Perform an explicit euler step to get guess
                // updates Y_i and N_i
                //this->explicitEuler(Y_i, N_i, dt);
                
                // N.M. equations value storage
                scalar fY(0.0);
                scalar fN(0.0);

                // Jacobian storage
                scalarSquareMatrix J(this->nEqns());
                scalarSquareMatrix invJ(this->nEqns());
                

                // temporary variables for Newton's method iterations
                scalar Y_next(0.0);
                scalar N_next(0.0);
                label converged(0);

                label newton_iteration(0);

                // Newton method iterations
                while( !converged && (newton_iteration <= 10)) 
                {
                    // get jacobian and inverse
                    this->jacobian(Y_i, N_i, dt, J, invJ);
                    
                    // newton method function values
                    this->newtonEquations(Y_i, N_i, Y_s[cellNumber], 
                    N_s[cellNumber], dt, fY, fN);

                    // Newton method step
                    Y_next = Y_i - (invJ(0,0) * fY + invJ(0,1) * fN);
                    N_next = N_i - (invJ(1,0) * fY + invJ(1,1) * fN);

                    Info << "Y_next, N_next: " << Y_next 
                        << ", " << N_next << endl;

                    // Force them to be zero if one is negative
                    if (Y_next <= 0. || N_next <= 0.)
                    {
                        Y_next = 0.0;
                        N_next = 0.0;
                    }
                    // else check for convergence
                    else if 
                        (fabs(Y_next - Y_i)/(Y_next + Foam::SMALL) < relTol 
                        && 
                        fabs(N_next - N_i)/(N_next + Foam::SMALL) < relTol)
                    {
                        converged = 1;
                    }

                    Y_i = Y_next;
                    N_i = N_next;

                    // increment the iteration counter
                    newton_iteration++;
                }// end Newton method loop

                // If we did not converge in 10 iterations
                if (newton_iteration > 10)
                {
                    // increment failure counter
                    newton_method_failures ++;

                    // FatalErrorInFunction
                    //     << "Newton Method did not converge" 
                    //         << abort(FatalError);

                    // if not converged set it as it was previously
                    // bad idea? 
                    Y_i = Y_0;
                    N_i = N_0;

                }
                
                // set source over the timestep for this cell
                this->Y_source[cellNumber] = -(Y_i - Y_s[cellNumber])/dt;
                this->N_source[cellNumber] = -(N_i - N_s[cellNumber])/dt;

            }// end for loop through cells
            
            Info << "The soot Newton's Method failed to converge in "
                << newton_method_failures << " cells." << endl;
        }
public: 

    // Constructor
    EulerImplicitSystem
    (
        const psiReactionThermo& thermo, 
        const basicSpecieMixture& composition,
        const volScalarField& Ns,
        const volScalarField& rho
    )
        :
        Y_source_dims(1,0,-1,0,0),// Source dimensions (exclude volume)
        N_source_dims(0,0,-1,0,0),// Source dimensions (exclude volume)
        MW_soot(composition.W(composition.species()["SOOT"])),
        MW_c2h2(composition.W(composition.species()["C2H2"])),
        MW_oh(composition.W(composition.species()["OH"])),
        MW_o2(composition.W(composition.species()["O2"])),
        rho(rho.primitiveField()),
        T(thermo.T().primitiveField()),
        Y_c2h2(composition.Y("C2H2").primitiveField()),
        Y_oh(composition.Y("OH").primitiveField()),
        Y_o2(composition.Y("O2").primitiveField()),
    Y_s(composition.Y("SOOT").primitiveField()), 
    N_s(Ns.primitiveField()),
    Y_source(rho.size()),
    N_source(rho.size()),
    C_A(0.0),
    C_1_1(0.0),
    C_1_2(0.0),
    C_1_3(0.0),
    C_1_4(0.0),
    C_2_1(0.0),
    C_2_2(0.0)
        {
            Info<< "Creating Soot Model Solver \n\n" << endl;
        }

    // Return the source for the Ysoot equation
    tmp<fvScalarMatrix> sourceY(const volScalarField& Y_soot)
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Y_soot, Y_source_dims));
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = this->Y_source;
            
            return tfvm;
        }
    
    // Return the source for the N_s equation
    tmp<fvScalarMatrix> sourceN(const volScalarField& N_soot)
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(N_soot, N_source_dims));

            // get reference to the fvmatrix inside the tmp
            fvScalarMatrix& fvm = tfvm.ref();

            fvm.source() = this->N_source;

            return tfvm;
        }

    // Main function for calculating the sources.
    // Manages substeps of the main time step dt
    void updateSources
    (
        const scalar dt,
        const label subSteps = 3, 
        const scalar relTol = 1e-6
    )
        {
            this->correct(dt, relTol);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
