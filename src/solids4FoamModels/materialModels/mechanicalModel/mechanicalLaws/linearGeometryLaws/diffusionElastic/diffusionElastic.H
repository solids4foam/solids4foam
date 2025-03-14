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
    diffusionElastic

Description
    Elastic mechanical law with the stress is given by

        sigma = mu*d*dev(symm(gradD)) + kappa*d*tr(gradD)*I

    where gradD is the displacement gradient, mu and kappa are user-defined
    scalar stiffness parameters and d is a diffusivity field specified by
    the "diffusivity" user setting, e.g. inverse distance from the wall.

    This law is useful as a mesh motion solver and also for checking the
    Jacobian in Newton procedures.

SourceFiles
    diffusionElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef diffusionElastic_H
#define diffusionElastic_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "motionDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class diffusionElastic Declaration
\*---------------------------------------------------------------------------*/

class diffusionElastic
:
    public mechanicalLaw
{
    // Private data

        //- Shear modulus type stiffness parameter
        const dimensionedScalar mu_;

        //- Bulk modulus type stiffness parameter
        const dimensionedScalar kappa_;

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
        diffusionElastic(const diffusionElastic&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionElastic&);

public:

    //- Runtime type information
    TypeName("diffusionElastic");

    // Constructors

        //- Construct from dictionary
        diffusionElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~diffusionElastic()
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
