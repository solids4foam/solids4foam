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
    Foam::fv::neoHookeanManufacturedSolutionSource

Description
    Applies the neo-Hookean manufactured-solution body force, per unit
    reference volume, to the momentum equation of the total Lagrangian solid
    model. OpenFOAM.com only; the other versions use
    neoHookeanManufacturedSolutionSolid.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    momentumSource
    {
        type             neoHookeanManufacturedSolutionSource;
        fields           (D);
    }
    \endverbatim

    The solution parameters are read from
    constant/neoHookeanManufacturedSolutionProperties.

Author
    Ivan Batistic, UCD.
    Pablo Castrillo, UCD.
    Philip Cardiff, UCD.

SourceFiles
    neoHookeanManufacturedSolutionSource.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanManufacturedSolutionSource_H
#define neoHookeanManufacturedSolutionSource_H

#ifdef OPENFOAM_COM

#include "fvOption.H"
#include "neoHookeanManufacturedSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class neoHookeanManufacturedSolutionSource Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanManufacturedSolutionSource
:
    public fv::option
{
    // Private Data

        //- Manufactured solution object
        neoHookeanManufacturedSolution mms_;


public:

    //- Runtime type information
    TypeName("neoHookeanManufacturedSolutionSource");


    // Constructors

        //- Construct from explicit source name and mesh
        neoHookeanManufacturedSolutionSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        neoHookeanManufacturedSolutionSource
        (
            const neoHookeanManufacturedSolutionSource&
        ) = delete;

        //- No copy assignment
        void operator=(const neoHookeanManufacturedSolutionSource&) = delete;


    //- Destructor
    virtual ~neoHookeanManufacturedSolutionSource() = default;


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

#endif

// ************************************************************************* //
