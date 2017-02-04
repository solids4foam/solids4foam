/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    neoHookeanElasticMisesPlasticRubin

Description
    Hyper-elasto-plastic mechanical law as described by Simo & Hughes (1998)
    in Box 9.1.

    Neo-Hookean elasticity with Mises/J2 plasticity.

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    More details found in:

    Simo & Hughes, Computational Inelasticity, 1998, Springer.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

SourceFiles
    neoHookeanElasticMisesPlasticRubin.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElasticMisesPlasticRubin_H
#define neoHookeanElasticMisesPlasticRubin_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElasticMisesPlasticRubin
:
    public mechanicalLaw
{
    // Private data

        //- Initial density
        const dimensionedScalar rho_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar k_;

        //- Kirchhoff yield stress
        volScalarField kappa_;

        //- Hardening modulus: do not confuse with the bulk modulus k
        dimensionedScalar K_;

        //- Relative deformation gradient volField
        volTensorField* relFPtr_;

        //- Jacobian of the toal deformation gradient volField
        volScalarField* JPtr_;

        //- Total deformation gradient volField
        volTensorField* FPtr_;

        //- Hydrostatic pressure == -hydrostatic stress
        volScalarField P_;

        //- Incremental change in sigmaY
        volScalarField Je_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        volSymmTensorField bEbarTrial_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        volSymmTensorField bEbar_;

        //- Plastic multiplier for Rubin algorithm
        volScalarField lambda_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- Switch between PC (plastic compaction) and C (J == Je) algorithms
        const Switch plasticCompaction_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        neoHookeanElasticMisesPlasticRubin
        (
            const neoHookeanElasticMisesPlasticRubin&
        );

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanElasticMisesPlasticRubin&);

        //- Make the relF field
        void makeRelF();

        //- Return a reference to the relF field
        volTensorField& relF();

        //- Make the relFf field
        void makeRelFf();

        //- Return a reference to the relFf field
        surfaceTensorField& relFf();

        //- Make the J field
        void makeJ();

        //- Return a reference to the J field
        volScalarField& J();

        //- Make the Jf field
        void makeJf();

        //- Return a reference to the Jf field
        surfaceScalarField& Jf();

        //- Make the F field
        void makeF();

        //- Return a reference to the F field
        volTensorField& F();

        //- Make the Ff field
        void makeFf();

        //- Return a reference to the Ff field
        surfaceTensorField& Ff();

        //- Calcualte Ibar such that det(bEbar) == 1
        tmp<volScalarField> Ibar
        (
            const volSymmTensorField& devBEbar
        );

public:

    //- Runtime type information
    TypeName("neoHookeanElasticMisesPlasticRubin");

    // Static data members


    // Constructors

        //- Construct from dictionary
        neoHookeanElasticMisesPlasticRubin
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~neoHookeanElasticMisesPlasticRubin();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();

        //- Return the desired new time-step
        virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
