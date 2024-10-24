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
    coconutCouplingInterface

Description
    This class ports the coupling procedures from the CoCoNuT toolbox available
    at https://github.com/pyfsi/coconut.


Author
    Ported to OpenFOAM by Philip Cardiff, UCD.
    See the CoCoNuT toolbox to see the authors or the original Python procedures.

SourceFiles
    coconutCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef coconutCouplingInterface_H
#define coconutCouplingInterface_H

#include "fluidSolidInterface.H"
#include "coconutModel.H"
#include "coconutConvergenceCriteria.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class coconutCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class coconutCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- Residual
        scalar residual_;

        //- Relaxation factor
        scalar relaxationFactor_;

        //- Coconuts model, e.g. least squares
        coconutModel model_;

        //- Coconuts convergence check
        coconutConvergenceCriteria convergenceCriteria_;

        //- Current interface values
        vectorField x_;

        //- Interface values from previous time steps
        PtrList<vectorField> xPrev_;


    // Private Member Functions

        //- Method to compute the residual (based on iqni.py logic)
        scalar computeResidual();

        //- Predictor
        tmp<vectorField> predictorPredict(const vectorField& x) const;

        //- Solver 0 solution step
        tmp<vectorField> solver0SolveSolutionStep(const vectorField& x);

        //- Solver 1 solution step
        tmp<vectorField> solver1SolveSolutionStep(const vectorField& x);

        //- Finialise the iteration
        void finalizeIteration(const vectorField& r);

        //- Disallow default bitwise copy construct
        coconutCouplingInterface(const coconutCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const coconutCouplingInterface&);


public:

    //- Runtime type information
    TypeName("coconutIQNILS");


    // Constructors

        //- Construct from components
        coconutCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~coconutCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Method to initialize the solver to be implemented in derived
            //  classes
            virtual void initialize();

            //- Method to perform one time step of the solver
            virtual void solveTimeStep();

            //- Method to restart from a given time step
            //virtual void restart(int timeStep);

            //- Method to finalize the solver
            virtual void finalize();

            //- Evolve the interface
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
