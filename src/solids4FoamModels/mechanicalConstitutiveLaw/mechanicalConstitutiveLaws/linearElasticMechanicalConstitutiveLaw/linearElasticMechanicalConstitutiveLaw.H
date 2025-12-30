/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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

Class
    Foam::linearElasticMechanicalConstitutiveLaw

Description
    Small-strain linear elastic mechanical constitutive law.

    This class implements an isotropic linear elastic constitutive
    relation in the small-strain regime, computing the Cauchy stress
    as a linear function of the symmetric displacement gradient.

    The law is stateless and does not require any constitutive history
    variables. All material parameters (e.g. Young's modulus and
    Poisson's ratio) are read from the supplied dictionary at
    construction time.

    Stress evaluation is performed using the
    smallStrainMechanicalConstitutiveLawKinematics interface. Finite-
    strain formulations are not supported by this law.

    Usage:
    The linear elastic law is selected via the mechanicalConstitutiveLaw
    run-time selection mechanism and used by the
    mechanicalConstitutiveLawManager during stress updates.

    Notes:
    - This law assumes infinitesimal strains and rotations.
    - No constitutive state variables are required.
    - Tangent-related quantities are not computed at this stage, but
      may be added later if required.


SourceFiles
    linearElasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMechanicalConstitutiveLaw_H
#define linearElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class linearElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private Data

        //- Young's modulus
        const dimensionedScalar E_;

        //- Poisson's ration
        const dimensionedScalar nu_;

        //- First Lamé parameter
        dimensionedScalar lambda_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;


    // Private Member Functions

        //- No copy construct
        linearElasticMechanicalConstitutiveLaw
        (
            const linearElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElastic");

    // Constructors

        //- Construct from dictionary
        explicit linearElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearElasticMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Edit

        //- Small-strain constitutive evaluation
        virtual void evaluate
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
