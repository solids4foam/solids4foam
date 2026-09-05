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

\*---------------------------------------------------------------------------*/

#include "leastSquaresReconstruction.H"
#include "compatibilityFunctions.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(leastSquaresReconstruction, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

leastSquaresReconstruction::leastSquaresReconstruction(const fvMesh& mesh)
:
#ifdef FOAMEXTEND
    MeshObject<fvMesh, leastSquaresReconstruction>(mesh),
#else
    MeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        leastSquaresReconstruction
    >(mesh),
#endif
    displacementSchemePtr_(),
    pressureSchemePtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

leastSquaresReconstruction::~leastSquaresReconstruction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

const leastSquaresScheme& leastSquaresReconstruction::scheme
(
    const word& fieldRole,
    const boolList& includePatchInStencils,
    const dictionary& dict
) const
{
    autoPtr<leastSquaresScheme>* schemePtr = nullptr;

    if (fieldRole == "displacement")
    {
        schemePtr = &displacementSchemePtr_;
    }
    else if (fieldRole == "pressure")
    {
        schemePtr = &pressureSchemePtr_;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown least-squares reconstruction field role '"
            << fieldRole << "'" << nl << nl
            << "Valid field roles are displacement and pressure"
            << exit(FatalError);
    }

    if (schemePtr->empty())
    {
        schemePtr->set
        (
            leastSquaresScheme::New
            (
                mesh(),
                includePatchInStencils,
                dict
            ).ptr()
        );
    }

    return autoPtrRef(*schemePtr);
}


void leastSquaresReconstruction::clear() const
{
    if (displacementSchemePtr_.valid())
    {
        displacementSchemePtr_->clear();
    }

    if (pressureSchemePtr_.valid())
    {
        pressureSchemePtr_->clear();
    }
}


#ifdef FOAMEXTEND
bool leastSquaresReconstruction::movePoints() const
#else
bool leastSquaresReconstruction::movePoints()
#endif
{
    clear();

    return true;
}


#ifdef FOAMEXTEND
bool leastSquaresReconstruction::updateMesh(const mapPolyMesh&) const
#else
void leastSquaresReconstruction::updateMesh(const mapPolyMesh&)
#endif
{
    displacementSchemePtr_.clear();
    pressureSchemePtr_.clear();

#ifdef FOAMEXTEND
    return true;
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
