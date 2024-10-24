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
    coconutConvergenceCriteria

Description
    This class ports the models from the CoCoNuT toolbox available
    at https://github.com/pyfsi/coconut.

Author
    Ported to OpenFOAM by Philip Cardiff, UCD.
    See the CoCoNuT toolbox to see the authors or the original Python procedures.

SourceFiles
    coconutConvergenceCriteria.C

\*---------------------------------------------------------------------------*/

#ifndef coconutConvergenceCriteria_H
#define coconutConvergenceCriteria_H

#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                        Class coconutConvergenceCriteria Declaration
\*---------------------------------------------------------------------------*/

class coconutConvergenceCriteria
{
    // Private data

        //- Tolerance
        const scalar tolerance_;

        //- Iteration
        int iter_;

        //- Maximum number of allowed iterations
        const int maxIter_;

        //- Initial norm
        scalar initialNorm_;

        //- Last norm
        scalar lastNorm_;

        //- Has the initial norm been set?
        bool isInitialNormSet_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coconutConvergenceCriteria(const coconutConvergenceCriteria&);

        //- Disallow default bitwise assignment
        void operator=(const coconutConvergenceCriteria&);


public:

    //- Runtime type information
    TypeName("coconutConvergenceCriteria");


    // Constructors

        //- Construct from components
        coconutConvergenceCriteria
        (
            const scalar tolerance,
            const int maxIter
        );


    // Destructor

        virtual ~coconutConvergenceCriteria()
        {}


    // Member Functions

        // Edit

            //- Initialise values at the start of a time step
            void initializeSolutionStep();

            //- Finalise values at the start of a time step
            void finaliseSolutionStep()
            {}

            //- Initialise values at the start of a time step
            void update(const vectorField& r);

            //- Initialise values at the start of a time step
            bool isSatisfied() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
