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

        // Pointer to the fluid mesh motion solver
        // We will use a solid solver for this
        mutable autoPtr<solidModel> meshMotionPtr_;

        //- Scale factor for the solid equations to preserve the condition
        //  number of the monlithic system
        //  We default this to be the reciprocal of the solid's shear modulus
        const scalar solidSystemScaleFactor_;

        //- Flag to enable fluid-to-solid coupling
        const Switch fluidToSolidCoupling_;

        //- Flag to enable solid-to-fluid coupling
        const Switch solidToFluidCoupling_;

        //- Number of sub-regions: 2 for fluid+solid, 3 for
        //  fluid+meshMotion+solid
        const label nRegions_;

        //- Pointer to sub-matrices
        Mat* subMatsPtr_;


    // Private Member Functions

        //- Make the fluid mesh motion solver
        void makeMeshMotion() const;

        //- Return the fluid mesh motion solver
        solidModel& meshMotion();

        //- Update mesh motion interface displacement
        //  This maps the current solid displacement interface field to the
        //  mesh motion displacement field. Both D and pointD are mapped
        void updateMeshMotionInterfaceDisplacement();

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
        {
            // Destroy each submatrix
            for (label i = 0; i < nRegions_*nRegions_; ++i)
            {
                MatDestroy(&subMatsPtr_[i]);
            }

            // Free the array of Mat handles
            delete[] subMatsPtr_;
            subMatsPtr_ = nullptr;
        }


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            virtual label initialiseJacobian(Mat& jac);

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x);

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
