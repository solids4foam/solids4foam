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
    Foam::fv::manufacturedSolutionSource

Group
    grpFvOptionsSources

Description
    Applies a manufactured solution vector body force to the momentum
    equation for solid solvers.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    manufacturedSolutionSource1
    {
        // Mandatory entries (unmodifiable)
        type             manufacturedSolutionSource;
        fields           (D);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

SourceFiles
    manufacturedSolutionSource.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef manufacturedSolutionSource_H
#define manufacturedSolutionSource_H

#include "fvOption.H"
#include "manufacturedSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class manufacturedSolutionSource Declaration
\*---------------------------------------------------------------------------*/

class manufacturedSolutionSource
:
    public fv::option
{
    // Private Data

        //- Manufactured solution object
        manufacturedSolution mms_;


public:

    //- Runtime type information
    TypeName("manufacturedSolutionSource");


    // Constructors

        //- Construct from explicit source name and mesh
        manufacturedSolutionSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        manufacturedSolutionSource(const manufacturedSolutionSource&) = delete;

        //- No copy assignment
        void operator=(const manufacturedSolutionSource&) = delete;


    //- Destructor
    virtual ~manufacturedSolutionSource() = default;


    // Member Functions

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            const volScalarField&,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
