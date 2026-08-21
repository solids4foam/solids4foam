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

#include "leastSquaresScheme.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(leastSquaresScheme, 0);
defineRunTimeSelectionTable(leastSquaresScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void leastSquaresScheme::makeFaceCentreValueCoeffs() const
{
    if
    (
        ownerFaceCentreValueCoeffsPtr_.valid()
     || neighbourFaceCentreValueCoeffsPtr_.valid()
    )
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    const CompactListList<label>& cellStencils = stencil().cellsStencil();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField& faceCentres = mesh_.Cf();

    labelList ownerRowSizes(mesh_.nFaces(), 0);
    labelList neighbourRowSizes(mesh_.nFaces(), 0);

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
    {
        ownerRowSizes[faceI] = cellStencils[owner[faceI]].size() + 1;
        neighbourRowSizes[faceI] =
            cellStencils[neighbour[faceI]].size() + 1;
    }

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        const labelUList& faceCells = mesh_.boundary()[patchI].faceCells();

        forAll(patch, patchFaceI)
        {
            ownerRowSizes[patch.start() + patchFaceI] =
                cellStencils[faceCells[patchFaceI]].size() + 1;
        }
    }

    ownerFaceCentreValueCoeffsPtr_.set
    (
        new CompactListList<scalar>(ownerRowSizes)
    );
    neighbourFaceCentreValueCoeffsPtr_.set
    (
        new CompactListList<scalar>(neighbourRowSizes)
    );

    CompactListList<scalar>& ownerCoeffs =
        *ownerFaceCentreValueCoeffsPtr_;
    CompactListList<scalar>& neighbourCoeffs =
        *neighbourFaceCentreValueCoeffsPtr_;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
    {
        UList<scalar> faceOwnerCoeffs = ownerCoeffs[faceI];
        cellValueCoeffsAtPoint
        (
            owner[faceI],
            faceCentres[faceI],
            faceOwnerCoeffs
        );

        UList<scalar> faceNeighbourCoeffs = neighbourCoeffs[faceI];
        cellValueCoeffsAtPoint
        (
            neighbour[faceI],
            faceCentres[faceI],
            faceNeighbourCoeffs
        );
    }

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        const labelUList& faceCells = mesh_.boundary()[patchI].faceCells();

        forAll(patch, patchFaceI)
        {
            const label faceI = patch.start() + patchFaceI;
            UList<scalar> faceOwnerCoeffs = ownerCoeffs[faceI];

            cellValueCoeffsAtPoint
            (
                faceCells[patchFaceI],
                faceCentres[faceI],
                faceOwnerCoeffs
            );
        }
    }
}


const CompactListList<scalar>&
leastSquaresScheme::ownerFaceCentreValueCoeffs() const
{
    if (ownerFaceCentreValueCoeffsPtr_.empty())
    {
        makeFaceCentreValueCoeffs();
    }

    return *ownerFaceCentreValueCoeffsPtr_;
}


const CompactListList<scalar>&
leastSquaresScheme::neighbourFaceCentreValueCoeffs() const
{
    if (neighbourFaceCentreValueCoeffsPtr_.empty())
    {
        makeFaceCentreValueCoeffs();
    }

    return *neighbourFaceCentreValueCoeffsPtr_;
}


template<class Type, class GlobalCells>
inline Type leastSquaresScheme::fieldValue
(
    const label globalCellID,
    const GlobalCells& globalCells,
    const Map<FixedList<label, 2>>& remoteLocation,
    const UList<Type>& localField,
    const List<Field<Type>>& remoteField
) const
{
    if (globalCells.isLocal(globalCellID))
    {
        return localField[globalCells.toLocal(globalCellID)];
    }

    if (debug && !remoteLocation.found(globalCellID))
    {
        FatalErrorInFunction
            << "Remote location mapping not found for global cell ID "
            << globalCellID << abort(FatalError);
    }

    const FixedList<label, 2>& loc = remoteLocation[globalCellID];
    return remoteField[loc[0]][loc[1]];
}


template<class Type, class GlobalCells>
Type leastSquaresScheme::evaluateCellValueCoeffs
(
    const label cellI,
    const UList<scalar>& coeffs,
    const GlobalCells& globalCells,
    const Map<FixedList<label, 2>>& remoteLocation,
    const UList<Type>& localField,
    const List<Field<Type>>& remoteField
) const
{
    const UList<label>& cellStencil = stencil().cellsStencil()[cellI];

    if (coeffs.size() != cellStencil.size() + 1)
    {
        FatalErrorInFunction
            << "Coefficient list size " << coeffs.size()
            << " does not match stencil plus owner size "
            << cellStencil.size() + 1 << " for cell " << cellI
            << abort(FatalError);
    }

    Type value = coeffs[cellStencil.size()]*localField[cellI];

    forAll(cellStencil, stencilI)
    {
        value +=
            coeffs[stencilI]
           *fieldValue
            (
                cellStencil[stencilI],
                globalCells,
                remoteLocation,
                localField,
                remoteField
            );
    }

    return value;
}


template<class Type>
void leastSquaresScheme::exchangeProcessorPatchField
(
    const processorFvPatch& procPatch,
    const Field<Type>& sendField,
    Field<Type>& receiveField
) const
{
    if (Pstream::myProcNo() < procPatch.neighbProcNo())
    {
        procPatch.send(Pstream::commsTypes::blocking, sendField);
        procPatch.receive(Pstream::commsTypes::blocking, receiveField);
    }
    else
    {
        procPatch.receive(Pstream::commsTypes::blocking, receiveField);
        procPatch.send(Pstream::commsTypes::blocking, sendField);
    }
}


template<class Type>
void leastSquaresScheme::evaluateFaceCentreValues
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& ownerValues,
    GeometricField<Type, fvsPatchField, surfaceMesh>& neighbourValues
) const
{
    const CompactListList<scalar>& ownerCoeffs =
        ownerFaceCentreValueCoeffs();
    const CompactListList<scalar>& neighbourCoeffs =
        neighbourFaceCentreValueCoeffs();
    const leastSquaresStencil& cellStencil = stencil();
    const globalIndex& globalCells = cellStencil.globalCells();
    const Map<FixedList<label, 2>>& remoteLocation =
        cellStencil.remoteCellLocation();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const UList<Type>& vfI = vf.internalField();
    const List<Field<Type>> remoteField =
        cellStencil.remoteFieldPerProc(vfI);
    Field<Type>& ownerValuesI = Foam::primitiveFieldRef(ownerValues);
    Field<Type>& neighbourValuesI = Foam::primitiveFieldRef(neighbourValues);

    if
    (
        ownerValuesI.size() != mesh_.nInternalFaces()
     || neighbourValuesI.size() != mesh_.nInternalFaces()
    )
    {
        FatalErrorInFunction
            << "Face-centre output fields do not match the mesh"
            << abort(FatalError);
    }

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
    {
        ownerValuesI[faceI] = evaluateCellValueCoeffs
        (
            owner[faceI],
            ownerCoeffs[faceI],
            globalCells,
            remoteLocation,
            vfI,
            remoteField
        );
        neighbourValuesI[faceI] = evaluateCellValueCoeffs
        (
            neighbour[faceI],
            neighbourCoeffs[faceI],
            globalCells,
            remoteLocation,
            vfI,
            remoteField
        );
    }

    auto& ownerBoundary = Foam::boundaryFieldRef(ownerValues);
    auto& neighbourBoundary = Foam::boundaryFieldRef(neighbourValues);

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const labelUList& faceCells = patch.faceCells();
        Field<Type>& patchOwnerValues = ownerBoundary[patchI];
        Field<Type>& patchNeighbourValues = neighbourBoundary[patchI];

        forAll(patch, patchFaceI)
        {
            const label faceI = patch.start() + patchFaceI;

            patchOwnerValues[patchFaceI] = evaluateCellValueCoeffs
            (
                faceCells[patchFaceI],
                ownerCoeffs[faceI],
                globalCells,
                remoteLocation,
                vfI,
                remoteField
            );
        }

        if (isA<processorFvPatch>(patch))
        {
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(patch);

            exchangeProcessorPatchField
            (
                procPatch,
                patchOwnerValues,
                patchNeighbourValues
            );
        }
        else
        {
            // There is no second cell reconstruction on a physical boundary.
            // Use a neutral value; boundary-condition policy belongs to the
            // caller, for example alphaStab.
            patchNeighbourValues = patchOwnerValues;
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<leastSquaresScheme> leastSquaresScheme::New
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
{
    const word schemeType
    (
        dict.lookupOrDefault<word>("type", "movingLeastSquares")
    );

    Info<< "Selecting least-squares reconstruction scheme "
        << schemeType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(schemeType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "leastSquaresScheme",
            schemeType,
            *dictionaryConstructorTablePtr_
        )   << exit(FatalIOError);
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(schemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("leastSquaresScheme::New(...)")
            << "Unknown least-squares reconstruction scheme type "
            << schemeType << nl << nl
            << "Valid least-squares reconstruction scheme types are:"
            << nl << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<leastSquaresScheme>
    (
        ctorPtr(mesh, includePatchInStencils, dict)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

leastSquaresScheme::leastSquaresScheme(const fvMesh& mesh)
:
    mesh_(mesh),
    ownerFaceCentreValueCoeffsPtr_(),
    neighbourFaceCentreValueCoeffsPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

leastSquaresScheme::~leastSquaresScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void leastSquaresScheme::clearFaceCentreValueCoeffs() const
{
    ownerFaceCentreValueCoeffsPtr_.clear();
    neighbourFaceCentreValueCoeffsPtr_.clear();
}


void leastSquaresScheme::faceCentreValues
(
    const volScalarField& vf,
    surfaceScalarField& ownerValues,
    surfaceScalarField& neighbourValues
) const
{
    evaluateFaceCentreValues<scalar>(vf, ownerValues, neighbourValues);
}


void leastSquaresScheme::faceCentreValues
(
    const volVectorField& vf,
    surfaceVectorField& ownerValues,
    surfaceVectorField& neighbourValues
) const
{
    evaluateFaceCentreValues<vector>(vf, ownerValues, neighbourValues);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
