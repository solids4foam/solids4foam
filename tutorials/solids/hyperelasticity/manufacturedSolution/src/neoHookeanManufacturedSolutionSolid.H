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
    Foam::solidModels::neoHookeanManufacturedSolutionSolid

Description
    Total Lagrangian solid with the neo-Hookean manufactured body force for
    OpenFOAM.org and foam-extend, where fvOptions are unavailable. The
    solution parameters are read from
    constant/neoHookeanManufacturedSolutionProperties, as for the OpenFOAM.com
    fvOption.

SourceFiles
    neoHookeanManufacturedSolutionSolid.C

\*----------------------------------------------------------------------------*/

#ifndef neoHookeanManufacturedSolutionSolid_H
#define neoHookeanManufacturedSolutionSolid_H

#ifndef OPENFOAM_COM

#include "nonLinGeomTotalLagTotalDispSolid.H"
#include "neoHookeanManufacturedSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                      Class neoHookeanManufacturedSolutionSolid Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanManufacturedSolutionSolid
:
    public nonLinGeomTotalLagTotalDispSolid
{
    // Private Data

        //- Manufactured solution supplying the body force
        neoHookeanManufacturedSolution mms_;


protected:

    //- Return the manufactured body force per unit reference volume
    virtual tmp<volVectorField> fvOptionsSource() const;


public:

    //- Runtime type information
    TypeName("neoHookeanManufacturedSolutionSolid");


    // Constructors

        //- Construct from time and region
        //  nonLinGeomTotalLagTotalDispSolid passes its own typeName to
        //  solidModel, so the coefficients are still read from
        //  nonLinearGeometryTotalLagrangianTotalDisplacementCoeffs
        neoHookeanManufacturedSolutionSolid(Time& runTime, const word& region);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif

// ************************************************************************* //
