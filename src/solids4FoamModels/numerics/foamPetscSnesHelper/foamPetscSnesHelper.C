/*---------------------------------------------------------------------------* \
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

#ifdef USE_PETSC

#include "foamPetscSnesHelper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(foamPetscSnesHelper, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamPetscSnesHelper::foamPetscSnesHelper
(
    const word& fieldName,
    fileName optionsFile,
    const fvMesh& mesh,
    const solutionLocation& location,
    const Switch stopOnPetscError,
    const Switch initialise,
    const fileName& fvSolutionLocation
)
:
    initialised_(initialise),
    fieldName_(fieldName),
    optionsFile_(optionsFile),
    fvSolutionLocation_(fvSolutionLocation),
    stopOnPetscError_(stopOnPetscError),
    assembly_(mesh, location),
    solverPtr_()
{
    if (initialise)
    {
        foamPetscSnesSolver::initialisePetsc();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

foamPetscSnesHelper::~foamPetscSnesHelper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

foamPetscSnesSolver& foamPetscSnesHelper::solver()
{
    if (solverPtr_.empty())
    {
        if (!initialised_)
        {
            FatalErrorInFunction
                << "The " << fieldName_ << " system cannot be solved because "
                << "the foamPetscSnesHelper object was not initialised during "
                << "construction" << abort(FatalError);
        }

        solverPtr_.set
        (
            new foamPetscSnesSolver
            (
                *this,
                fieldName_,
                optionsFile_,
                assembly_.mesh(),
                stopOnPetscError_,
                fvSolutionLocation_
            )
        );
    }

    return solverPtr_();
}


label foamPetscSnesHelper::initialiseJacobian
(
    Mat& jac,
    const fvMesh& mesh,
    const label blockSize,
    const bool createMat
) const
{
    if (&mesh != &assembly_.mesh())
    {
        FatalErrorInFunction
            << "The mesh " << mesh.name() << " is not the mesh of this "
            << "helper: use the helper of the mesh" << abort(FatalError);
    }

    return assembly_.initialiseJacobian(jac, blockSize, createMat);
}


label foamPetscSnesHelper::initialiseSolution
(
    Vec& x,
    const fvMesh& mesh,
    const label blockSize,
    const bool createVec
) const
{
    if (&mesh != &assembly_.mesh())
    {
        FatalErrorInFunction
            << "The mesh " << mesh.name() << " is not the mesh of this "
            << "helper: use the helper of the mesh" << abort(FatalError);
    }

    return assembly_.initialiseSolution(x, blockSize, createVec);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
