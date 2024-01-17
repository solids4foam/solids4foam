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

#include "dualMechanicalModel.H"
#include "fvc.H"
#include "fvcGradf.H"
#include "gaussGrad.H"
#include "twoDPointCorrector.H"
#include "fixedGradientFvPatchFields.H"
#include "wedgePolyPatch.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "ZoneIDs.H"
#else
    #include "ZoneID.H"
    #include "crackerFvMesh.H"
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dualMechanicalModel::makeCellInThisMaterialList() const
{
    if (!cellInThisMaterialList_.empty())
    {
        FatalError
            << "Pointers already set!" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;
    cellInThisMaterialList_.resize(laws.size());

    // Take a reference to the primary mesh
    const fvMesh& primaryMesh = mechModel_.mesh();

    // Create field with material law index, i.e. a cell in material 3 is
    // assigned the scalar 3
    volScalarField lawIDField
    (
        IOobject
        (
            "lawIDField",
            primaryMesh.time().timeName(),
            primaryMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh,
        dimensionedScalar("0", dimless, 0)
    );
    const scalarField& lawIDFieldI = lawIDField.internalField();

    PtrList<volScalarField> subMeshLawIDFields
    (
        mechModel_.solSubMeshes().subMeshes().size()
    );

    forAll(subMeshLawIDFields, subMeshI)
    {
        // Create lawID for this subMesh
        // Note that lawI and subMeshI are the same
        subMeshLawIDFields.set
        (
            subMeshI,
            new volScalarField
            (
                IOobject
                (
                    "subMeshLawIDField",
                    primaryMesh.time().timeName(),
                    mechModel_.solSubMeshes().subMeshes()[subMeshI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mechModel_.solSubMeshes().subMeshes()[subMeshI].subMesh(),
                dimensionedScalar("lawI", dimless, subMeshI),
                "zeroGradient"
            )
        );
    }

    // Populate global field with subMesh field
    mechModel_.solSubMeshes().mapSubMeshVolFields<scalar>
    (
        subMeshLawIDFields, lawIDField
    );

    forAll(laws, lawI)
    {
        cellInThisMaterialList_.set
        (
            lawI,
            new volScalarField
            (
                IOobject
                (
                    "cellInThisMaterial" + Foam::name(lawI),
                    primaryMesh.time().timeName(),
                    primaryMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                primaryMesh,
                dimensionedScalar("0", dimless, 0),
                "zeroGradient"
            )
        );

        scalarField& cellInThisMaterialListI = cellInThisMaterialList_[lawI];
        forAll(lawIDFieldI, cellI)
        {
            if (mag(lawIDFieldI[cellI] - lawI) < SMALL)
            {
                cellInThisMaterialListI[cellI] = 1;
            }
        }
        cellInThisMaterialList_[lawI].correctBoundaryConditions();
    }
}


const Foam::PtrList<Foam::volScalarField>&
Foam::dualMechanicalModel::cellInThisMaterialList() const
{
    if (cellInThisMaterialList_.empty())
    {
        makeCellInThisMaterialList();
    }

    return cellInThisMaterialList_;
}


void Foam::dualMechanicalModel::makeDualFaceInThisMaterialList() const
{
    if (!dualFaceInThisMaterialList_.empty())
    {
        FatalError
            << "Pointers already set!" << abort(FatalError);
    }

    const PtrList<mechanicalLaw>& laws = *this;
    dualFaceInThisMaterialList_.resize(laws.size());

    forAll(laws, lawI)
    {
        dualFaceInThisMaterialList_.set
        (
            lawI,
            new surfaceScalarField
            (
                IOobject
                (
                    "dualFaceInThisMaterial" + Foam::name(lawI),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("0", dimless, 0)
            )
        );

        // For all dual faces, find the corresponding primary mesh cell
        // then determine which material it belongs to using the
        // cellInThisMaterial field from the primary mechanicalModel
        scalarField& dualFaceMask = dualFaceInThisMaterialList_[lawI];
        const scalarField& cellMask =
            cellInThisMaterialList()[lawI].internalField();
        forAll(dualFaceMask, dualFaceI)
        {
            // Primary mesh cell
            const label cellID = dualFaceToCell_[dualFaceI];

            // Material in primary mesh cell
            const scalar cellInThisMat = cellMask[cellID];

            dualFaceMask[dualFaceI] = cellInThisMat;
        }

        /*
        // Boundary faces are not set, as the dualFaceToCell map may not be
        // defined. In any case, it shouldn't be needed
        forAll(dualFaceInThisMaterialList_[lawI].boundaryField(), patchI)
        {
            scalarField& dualFaceMaskP =
                dualFaceInThisMaterialList_[lawI].boundaryFieldRef()[patchI];

            forAll(dualFaceMaskP, dualFaceI)
            {
                // Dual face ID
                const label dFaceID =
                    mesh_.boundaryMesh()[patchI].start() + dualFaceI;

                // Primary mesh cell
                const label cellID = dualFaceToCell_[dFaceID];

                // Material in primary mesh cell
                const scalar cellInThisMat = cellMask[cellID];

                dualFaceMaskP[dualFaceI] = cellInThisMat;
            }
        }
        */
    }
}


const Foam::PtrList<Foam::surfaceScalarField>&
Foam::dualMechanicalModel::dualFaceInThisMaterialList() const
{
    if (dualFaceInThisMaterialList_.empty())
    {
        makeDualFaceInThisMaterialList();
    }

    return dualFaceInThisMaterialList_;
}


void Foam::dualMechanicalModel::clearOut()
{
    // Clear the list of mechanical laws
    PtrList<mechanicalLaw>::clear();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualMechanicalModel::dualMechanicalModel
(
    const fvMesh& mesh,
    const nonLinearGeometry::nonLinearType& nonLinGeom,
    const bool incremental,
    const mechanicalModel& mechModel,
    const labelList& dualFaceToCell
)
:
    PtrList<mechanicalLaw>(),
    mesh_(mesh),
    mechModel_(mechModel),
    planeStress_(mechModel.lookup("planeStress")),
    incremental_(incremental),
    dualFaceToCell_(dualFaceToCell),
    cellInThisMaterialList_(),
    dualFaceInThisMaterialList_()
{
    Info<< "Creating the dualMechanicalModel" << endl;

    // Read the mechanical laws
    PtrList<entry> lawEntries(mechModel.lookup("mechanical"));

    PtrList<mechanicalLaw>& laws = *this;
    laws.setSize(lawEntries.size());

    // Create mechanical laws
    forAll(laws, lawI)
    {
        if (nonLinGeom == nonLinearGeometry::LINEAR_GEOMETRY)
        {
            laws.set
            (
                lawI,
                mechanicalLaw::NewLinGeomMechLaw
                (
                    lawEntries[lawI].keyword(),
                    mesh,
                    lawEntries[lawI].dict(),
                    nonLinGeom
                )
            );
        }
        else if
        (
            nonLinGeom == nonLinearGeometry::UPDATED_LAGRANGIAN
         || nonLinGeom == nonLinearGeometry::TOTAL_LAGRANGIAN
        )
        {
            laws.set
            (
                lawI,
                mechanicalLaw::NewNonLinGeomMechLaw
                (
                    lawEntries[lawI].keyword(),
                    mesh,
                    lawEntries[lawI].dict(),
                    nonLinGeom
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "Foam::dualMechanicalModel::dualMechanicalModel(...)\n"
            )   << "It is not clear what type of mechanical law should be "
                << "created for a solidModel with nonLinGeom = " << nonLinGeom
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dualMechanicalModel::~dualMechanicalModel()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::dualMechanicalModel::mesh() const
{
    return mesh_;
}

#ifdef OPENFOAM_NOT_EXTEND
Foam::tmp<Foam::Field<Foam::RectangularMatrix<Foam::scalar>>>
Foam::dualMechanicalModel::materialTangentFaceField() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    // Prepare the field
    tmp< Field<RectangularMatrix<scalar>> > tresult
    (
        new Field<RectangularMatrix<scalar>>(mesh().nFaces(), RectangularMatrix<scalar>(6))
    );
    Field<RectangularMatrix<scalar>>& result = tresult.ref();

    if (laws.size() == 1)
    {
        result = laws[0].materialTangentField();

        if (result.size() != mesh().nFaces())
        {
            FatalErrorIn("dualMechanicalModel::materialTangentField()")
                << "The materialTangentField field for law 0 is the wrong size!"
                << abort(FatalError);
        }
    }
    else
    {
        // Accumulate data for all materials
        // Note: the value on each dual face is uniquely set by one material law
        // forAll(laws, lawI)
        // {
        //     const Field<scalarSquareMatrix> matTanI
        //     (
        //         laws[lawI].materialTangentField()
        //     );

        //     if (matTanI.size() != mesh().nFaces())
        //     {
        //         FatalErrorIn("dualMechanicalModel::materialTangentField()")
        //             << "The materialTangentField field for law " << lawI
        //             << " is the wrong size!" << abort(FatalError);
        //     }

        //     // Insert values from actual material region into main sigma field
        //     result += dualFaceInThisMaterialList()[lawI]*matTanI;
        // }

        // This does work but the problem can be that the maps for the boundary
        // faces sometimes struggle to be defined when creating the
        // dualFaceInThisMaterialList function. It should be possible to make
        // this robust
        notImplemented("Not implemented for more than one material");
    }

    return tresult;
}
#endif

Foam::tmp<Foam::surfaceScalarField>
Foam::dualMechanicalModel::bulkModulus() const
{
    const PtrList<mechanicalLaw>& laws = *this;

    // Prepare the field
    tmp< surfaceScalarField > tresult
    (
        new surfaceScalarField
        (
            IOobject
            (
                "dualBulkModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    );
#ifdef OPENFOAM_NOT_EXTEND
    surfaceScalarField& result = tresult.ref();
#else
    surfaceScalarField& result = tresult();
#endif

    if (laws.size() == 1)
    {
        result = fvc::interpolate(laws[0].bulkModulus());
    }
    else
    {
        // Accumulate data for all fields
        // Each face in the dualMesh lies in one cell (and hence one material)
        // in the primary mesh
        forAll(laws, lawI)
        {
            // Insert values from actual material region into main sigma field
            result +=
                dualFaceInThisMaterialList()[lawI]*fvc::interpolate
                (
                    laws[lawI].bulkModulus()
                );
        }
    }

    //result.write();

    return tresult;
}

void Foam::dualMechanicalModel::correct(surfaceSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    if (laws.size() == 1)
    {
        laws[0].correct(sigma);
    }
    else
    {
        // Accumulate data for all fields
        // Each face in the dualMesh lies in one cell (and hence one material)
        // in the primary mesh
        // The current approach is that each material calculates the stress
        // for all dual faces, and then the values of the faces that are
        // actually in that material are inserted into the global stress to be
        // returned to the solver

        // Reset stress to zero then accumulate it
        sigma = dimensionedSymmTensor("0", dimPressure, symmTensor::zero);

        forAll(laws, lawI)
        {
            // Create temporary stress field for this material
            surfaceSymmTensorField curSigma
            (
                IOobject
                (
                    "sigmafDualMeshLaw" + Foam::name(lawI),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            );

            // Calculate stress for all faces using this material law
            laws[lawI].correct(curSigma);

            // Insert values from actual material region into main sigma field
            sigma += dualFaceInThisMaterialList()[lawI]*curSigma;
        }
    }
}


// ************************************************************************* //
