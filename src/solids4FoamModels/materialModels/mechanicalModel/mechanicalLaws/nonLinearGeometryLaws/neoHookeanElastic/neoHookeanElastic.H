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
    neoHookeanElastic

Description
    Neo-Hookean hyper-elasticity, similar to model proposed in Simo & Hughes
    (1998) in Eqn 9.2.6, where the Cauchy stress is given as:

        sigma = 0.5*K*I*(J^2 - 1)/J + mu*dev(J^(-2.0/3.0)*B)/J

    where
    K        bulk modulus
    J        Jacobian of the deformation gradient
    I        indentity tensor
    mu       shear modulus
    B        left Cauchy–Green deformation tensor = F & F.T()
    F        deformation gradient

    If thermal stresses are active then 3.0*K*alpha*deltaT*I is added.

    Either E (Young's modulus) and nu (Poisson's ratio) can be provided or mu
    (shear modulus) and K (bulk modulus)

    Reference
    Simo & Hughes, Computational Inelasticity, 1998, Springer.


SourceFiles
    neoHookeanElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElastic_H
#define neoHookeanElastic_H

#include "mechanicalLaw.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "surfaceFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElastic
:
    public mechanicalLaw
{
    // Private data

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Alternative definition of the Kirchhoff pressure
        //  default: -0.5*K*(pow(J, 2.0) - 1.0)
        //  alternate: -0.5*K*(J - 1.0)
        const Switch alternatePressureDefinition_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        neoHookeanElastic(const neoHookeanElastic&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanElastic&);

public:

    //- Runtime type information
    TypeName("neoHookeanElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        neoHookeanElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~neoHookeanElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

#ifdef OPENFOAM_COM
        //- Return material tangent
        tmp<Field<scalarSquareMatrix>> materialTangentField() const;
#endif

        //- Return the bulk modulus field
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return the shear modulus field
        virtual tmp<volScalarField> shearModulus() const;

        //- Calculate the volField stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the surfaceField stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Calculate the stress using the displacement gradient
        virtual void correct
        (
            surfaceSymmTensorField& sigma,
            const surfaceTensorField& gradD
        ) const;

        //- Calculate the stress using the deformation gradient
        virtual void correctF
        (
            surfaceSymmTensorField& sigma,
            const surfaceTensorField& F
        ) const;

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
