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
    HolzapfelGasserOgdenElastic

Description
    Incompressible Holzapfel-Gasser-Ogden anisotropic hyperelastic law for
    coupledPressureDisplacementSolid pressure-displacement formulations.

    The law requires the following entries in mechanicalProperties:

    \verbatim
        type            HolzapfelGasserOgdenElastic;
        rho             rho [1 -3 0 0 0 0 0] <value>;
        E               E   [1 -1 -2 0 0 0 0] <value>;
        nu              nu  [0 0 0 0 0 0 0] 0.5;
        // or
        mu              mu  [1 -1 -2 0 0 0 0] <value>;

        fibreAngle      fibreAngle [0 0 0 0 0 0 0] <degrees>;
        k1              k1 [1 -1 -2 0 0 0 0] <value>;
        k2              k2 [0 0 0 0 0 0 0] <value>;
        impKcoeff       <optional scalar>;
    \endverbatim

    The local cylindrical material directions must be present as volVectorFields
    Ec, Ea and Er, and as surfaceVectorFields Ecf, Eaf and Erf. For cases such
    as arterial geometries these fields are typically generated before running
    the solver.

SourceFiles
    HolzapfelGasserOgdenElastic.C

Author
    Zeljko Tukovic, FSB Zagreb.

\*---------------------------------------------------------------------------*/

#ifndef HolzapfelGasserOgdenElastic_H
#define HolzapfelGasserOgdenElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class HolzapfelGasserOgdenElastic Declaration
\*---------------------------------------------------------------------------*/

class HolzapfelGasserOgdenElastic
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

        //- Calculate the stress on faces (shared by correct(surfaceSymmTensorField&))
        void correctF
        (
            surfaceSymmTensorField& sigma,
            const surfaceTensorField& Ff
        );

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
        HolzapfelGasserOgdenElastic(const HolzapfelGasserOgdenElastic&);

        //- Disallow default bitwise assignment
        void operator=(const HolzapfelGasserOgdenElastic&);

public:

    //- Runtime type information
    TypeName("HolzapfelGasserOgdenElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        HolzapfelGasserOgdenElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~HolzapfelGasserOgdenElastic();


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
