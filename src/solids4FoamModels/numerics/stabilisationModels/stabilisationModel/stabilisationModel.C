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

InClass
    Foam::stabilisationModel

\*---------------------------------------------------------------------------*/

#include "stabilisationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stabilisationModel, 0);
    defineRunTimeSelectionTable(stabilisationModel, stabModel);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stabilisationModel::stabilisationModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    mesh_(mesh),
    dict_(dict),
    dims_(dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
    scaleFactorJacobian_
    (
        dict.lookupOrDefault<scalar>("scaleFactorJacobian", scaleFactor_)
    ),
    faceScalarPtr_(),
    faceVectorPtr_(),
    scalarJacobianPtr_(),
    vectorJacobianPtr_(),
    h2Ptr_()
{}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::stabilisationModel> Foam::stabilisationModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
{
    const word modelType(dict.lookup("type"));

    Info<< "Selecting stabilisation model " << modelType << endl;

    auto* ctorPtr = stabModelConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "stabilisationModel",
            modelType,
            *stabModelConstructorTablePtr_
        )   << exit(FatalIOError);
    }

    return autoPtr<stabilisationModel>(ctorPtr(mesh, dict, dims));
}


// ************************************************************************* //
