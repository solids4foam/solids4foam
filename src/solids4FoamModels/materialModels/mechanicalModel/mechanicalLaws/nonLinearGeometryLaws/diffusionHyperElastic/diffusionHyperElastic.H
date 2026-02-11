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
    diffusionHyperElastic

Description
    Neo-Hookean hyper-elasticity where the Cauchy stress is given as:

        sigma = 0.5*K*I*(J^2 - 1)/J + mu*dev(J^(-2.0/3.0)*B)/J

    and mu and K are user-defined scalar stiffness parameters and d is a
    diffusivity field specified by the "diffusivity" user setting, e.g. inverse
    distance from the wall.

    This law is useful as a mesh motion solver.

SourceFiles
    diffusionHyperElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef diffusionHyperElastic_H
#define diffusionHyperElastic_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "motionDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class diffusionHyperElastic Declaration
\*---------------------------------------------------------------------------*/

class diffusionHyperElastic
:
    public mechanicalLaw
{
    // Private data

        //- Shear modulus type stiffness parameter
        const dimensionedScalar mu_;

        //- Bulk modulus type stiffness parameter
        const dimensionedScalar K_;

        //- Maximum factor
        //  The max value of d is limit to maxFactor*average(d)
        //  Defaults to 10
        const scalar maxFactor_;

        //- Minimum factor
        //  The min value of d is limit to minFactor*average(d)
        //  Defaults to 0.1
        const scalar minFactor_;

        //- Stiffness calculator, which calculates d
        mutable autoPtr<motionDiffusivity> motionDiffPtr_;

        //- d surface field
        mutable surfaceScalarField df_;

        //- d vol field
        mutable volScalarField d_;

    // Private Member Functions

        //- Update the d fields
        void updateD() const;

        //- Disallow default bitwise copy construct
        diffusionHyperElastic(const diffusionHyperElastic&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionHyperElastic&);

public:

    //- Runtime type information
    TypeName("diffusionHyperElastic");

    // Constructors

        //- Construct from dictionary
        diffusionHyperElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~diffusionHyperElastic()
        {}


    // Member Functions

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return shear modulus
        virtual tmp<volScalarField> shearModulus() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
