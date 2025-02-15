/*---------------------------------------------------------------------------*\
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
    GuccioneElastic

Description
    Guccione hyper-elasticity law suitable for myocardium. This implementation
    is based on the approach described in E Garcia-Blanco, R Ortigosa, AJ Gil,
    CH Lee, J Bonet, A new computational framework for electro-activation in
    cardiac mechanics, Comput. Methods Appl. Mech. Engrg. 348 (2019) 796–845.

SourceFiles
    GuccioneElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef GuccioneElastic_H
#define GuccioneElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class GuccioneElastic
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

        //- Flag to indicate a uniform fibre field
        //  If true, the fibre direction is read from mechanicalProperties
        //  If false, the fibre field is read from the starting time directory
        const Switch uniformFibreField_;

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
        surfaceSymmTensorField Sf_;


    // Private Member Functions

        //- Create f0 field
        tmp<volVectorField> makeF0
        (
            const Switch& uniformFibreField,
            const fvMesh& mesh,
            const dictionary& dict
        ) const;

        //- Create f0f field
        tmp<surfaceVectorField> makeF0f
        (
            const Switch& uniformFibreField,
            const fvMesh& mesh,
            const dictionary& dict
        ) const;

        //- Calculate stress given gradD and pressure
        void calculateStress
        (
            surfaceSymmTensorField& sigma,
            const surfaceTensorField& gradD
        );

        //- Disallow default bitwise copy construct
        GuccioneElastic(const GuccioneElastic&);

        //- Disallow default bitwise assignment
        void operator=(const GuccioneElastic&);

public:

    //- Runtime type information
    TypeName("GuccioneElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        GuccioneElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~GuccioneElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

#ifdef OPENFOAM_NOT_EXTEND
        //- Material tangent for Newton Raphson loop
        virtual tmp<Field<scalarSquareMatrix>> materialTangentField() const;
#endif

        //- Return the bulk modulus field
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return the shear modulus field
        virtual tmp<volScalarField> shearModulus() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
