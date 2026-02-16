/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    pdGuccioneElastic

Description
    Guccione hyper-elasticity law suitable for myocardium. This implementation
    is based on the approach described in E Garcia-Blanco, R Ortigosa, AJ Gil,
    CH Lee, J Bonet, A new computational framework for electro-activation in
    cardiac mechanics, Comput. Methods Appl. Mech. Engrg. 348 (2019) 796–845.

SourceFiles
    pdGuccioneElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB

\*---------------------------------------------------------------------------*/

#ifndef pdGuccioneElastic_H
#define pdGuccioneElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class pdGuccioneElastic
:
    public mechanicalLaw
{
    // Private data

        //- Bulk modulus
        dimensionedScalar bulkModulus_;

        //- Guccione modulus parameter
        dimensionedScalar k_;

        //- Material parameters
        const scalar cf_;
        const scalar ct_;
        const scalar cfs_;

        //- Shear modulus in the limited of small strains
        dimensionedScalar mu_;

        //- Effective shear modulus
        volScalarField muEff_;

        //- Fibre directions in the reference configuration
        volVectorField f0_;
        surfaceVectorField f0f_;

        //- Sheet directions in the reference configuration
        //  f0 & s0 == 0
        volVectorField s0_;
        surfaceVectorField s0f_;

        //- Sheet normal directions in the reference configuration
        //  n0 & f0 == n0 & s0 == 0
        volVectorField n0_;
        surfaceVectorField n0f_;

        //- Rotation matrix from local fibre coordinate system to the global
        //  coordinate system
        volTensorField R_;
        surfaceTensorField Rf_;

        //- Outer product of the fibre directions with themselves
        volSymmTensorField f0f0_;
        surfaceSymmTensorField f0f0f_;

        //- 2nd Piola-Kirchhoff stress (without the hydrostatic term)
        volSymmTensorField S_;

        //- 2nd Piola-Kirchhoff stress (without the hydrostatic term)
        surfaceSymmTensorField Sf_;

        //- Exponential of Qf
        surfaceScalarField expQf_;

        //- Scaling coeff for impK
        const scalar impKcoeff_;

    // Private Member Functions

        //- Calculate Cauchy stress as a function of deformation gradient
        void calcDevCauchy
        (
            const tensor& F,
            const symmTensor& f0f0,
            const tensor& R,
            symmTensor& devCauchy
        ) const;

        // Calculate initial shear modulus by F perturabation
        void calcInitShearModulus();

        // Calculate effective shear modulus
        void calcEffectiveShearModulus();

        //- Disallow default bitwise copy construct
        pdGuccioneElastic(const pdGuccioneElastic&);

        //- Disallow default bitwise assignment
        void operator=(const pdGuccioneElastic&);

public:

    //- Runtime type information
    TypeName("pdGuccioneElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        pdGuccioneElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~pdGuccioneElastic();


    // Member Functions

        //- Shear modulus
        virtual const dimensionedScalar& mu() const
        {
            return mu_;
        }

        virtual const dimensionedScalar& k() const
        {
            return k_;
        }

        virtual const scalar& cf() const
        {
            return cf_;
        }

        virtual const scalar& cfs() const
        {
            return cfs_;
        }

        virtual const scalar& ct() const
        {
            return ct_;
        }

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
