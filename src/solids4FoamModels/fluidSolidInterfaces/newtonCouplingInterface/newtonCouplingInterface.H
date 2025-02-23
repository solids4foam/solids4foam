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
    newtonCouplingInterface

Description
    Monolithic fluid-solid interaction solver with Newton-Raphson solution
    algorithm.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    newtonCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef newtonCouplingInterface_H
#define newtonCouplingInterface_H

#include "fluidSolidInterface.H"
#include "foamPetscSnesHelper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class newtonCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class newtonCouplingInterface
:
    public fluidSolidInterface,
    public foamPetscSnesHelper
{
    // Private data

        //- My data
        // const scalar relaxationFactor_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        newtonCouplingInterface(const newtonCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const newtonCouplingInterface&);


public:

    //- Runtime type information
    TypeName("NewtonMonolithic");


    // Constructors

        //- Construct from components
        newtonCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~newtonCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Initialise the Jacobian matrix
            //  This procedure should set the non-zero structure of the matrix
            virtual label initialiseJacobian(Mat jac);

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector
            virtual label formResidual
            (
                PetscScalar *f,        // Residual
                const PetscScalar *x   // Solution
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,              // Jacobian
                const PetscScalar *x  // Solution
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
