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
    pdHGO

Description

SourceFiles
    pdHGO.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef pdHGO_H
#define pdHGO_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class pdHGO
:
    public mechanicalLaw
{
    // Private data

        // Initial density
        const dimensionedScalar rho_;

        // Shear modulus
        dimensionedScalar mu_;

        //- Effective shear modulus
        volScalarField muEff_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Fibre angle
        dimensionedScalar fibreAngle_;

        //- HGO constant
        dimensionedScalar k1_;

        //- HGO constant
        dimensionedScalar k2_;

        //- Circumferential direction vector
        volVectorField Ec_;

        //- Axial direction vector
        volVectorField Ea_;

        //- Radial direction vector
        volVectorField Er_;

        //- Circumferential direction vector
        surfaceVectorField Ecf_;

        //- Axial direction vector
        surfaceVectorField Eaf_;

        //- Radial direction vector
        surfaceVectorField Erf_;

        //- Scaling coeff for impK
        const scalar impKcoeff_;

    // Private Member Functions

        //- Calculate Cauchy stress as a function of deformation gradient
        void calcDevCauchy
        (
            const tensor& F,
            const vector& Ec,
            const vector& Ea,
            const vector& Er,
            symmTensor& devCauchy
        ) const;

        // Calculate initial shear modulus by F perturabation
        void calcInitialShearModulus();

        // Calculate effective shear modulus
        void calcEffectiveShearModulus();

        //- Disallow default bitwise copy construct
        pdHGO(const pdHGO&);

        //- Disallow default bitwise assignment
        void operator=(const pdHGO&);

public:

    //- Runtime type information
    TypeName("pdHGO");

    // Static data members


    // Constructors

        //- Construct from dictionary
        pdHGO
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~pdHGO();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart();

        //- Update total accumulated fields
        virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
