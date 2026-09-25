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
    manufacturedSolution

Description
    Class that generates a manufactured displacement solution and related stress
    solution and source term for a solid solver.

    The manufactured solution for displacement takes the form:

    D =
        vector
        (
            ax*sin(4*pi*x)*sin(2*pi*y)*sin(pi*z),
            ay*sin(4*pi*x)*sin(2*pi*y)*sin(pi*z),
            az*sin(4*pi*x)*sin(2*pi*y)*sin(pi*z)
        )

    where ax, ay, and az are user provided scalar parameters.

    Linear elasticity is assumed for calculating the stress and momentum
    equation source term, hence the Young's modulus (E) and Poisson's
    ratio (nu) must be provided.

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

SourceFiles
    manufacturedSolution.C

\*---------------------------------------------------------------------------*/

#ifndef manufacturedSolution_H
#define manufacturedSolution_H

#include "point.H"
#include "symmTensor.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class manufacturedSolution Declaration
\*---------------------------------------------------------------------------*/

class manufacturedSolution
{
    // Private data

        //- Const reference to the mesh (need for the body forces field)
        const fvMesh& mesh_;

        //- User parameters
        const scalar ax_;
        const scalar ay_;
        const scalar az_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Lame parameters
        const scalar mu_;
        const scalar lambda_;

        //- Body forces
        mutable autoPtr<volVectorField> bodyForcesPtr_;

    // Private Member Functions

        //- Calculate the body forces field
        void calcBodyForces() const;

        //- Disallow default bitwise copy construct
        manufacturedSolution
        (
            const manufacturedSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const manufacturedSolution&);


public:

    //- Runtime type information
    TypeName("manufacturedSolution");


    //- Constructor

        //- Construct from mesh and parameters
        manufacturedSolution
        (
            const fvMesh& mesh,
            const scalar ax,
            const scalar ay,
            const scalar az,
            const scalar E,
            const scalar nu
        );

        //- Construct from mesh and dictionary
        manufacturedSolution
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~manufacturedSolution()
    {}


    // Member Functions

        //- Calculate the stress at the given point
        symmTensor calculateStress(const vector& point);

        //- Calculate the strain at the given point
        symmTensor calculateStrain(const vector& point);

        //- Calculate the stress at the given point
        vector calculateDisplacement(const vector& point);

        //- Body forces
        const volVectorField& bodyForces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
