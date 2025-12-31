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

#include "mechanicalConstitutiveLawManager.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::mechanicalConstitutiveLawManager
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    laws_(),
    states_(),
    boundaryStates_(),
    lawCells_(),
    lawBoundaryFaces_()
{
    // Read the mechanical laws
    const PtrList<entry> lawEntries(dict.lookup("mechanical"));

    laws_.setSize(lawEntries.size());
    states_.setSize(lawEntries.size());
    boundaryStates_.setSize(lawEntries.size());
    lawCells_.setSize(lawEntries.size());
    lawBoundaryFaces_.setSize(lawEntries.size());

    // Look list of law names
    wordList lawNames(lawEntries.size());
    forAll(lawNames, lawI)
    {
        lawNames[lawI] = lawEntries[lawI].keyword();
    }

    // Create a map for each cell to its mechanical law
    labelList cellToLaw(mesh_.nCells(), -1);

    forAll(lawNames, lawI)
    {
        const word& lawName = lawNames[lawI];
        const dictionary& lawDict = lawEntries[lawI].dict();

        // Construct law
        laws_.set
        (
            lawI,
            mechanicalConstitutiveLaw::New(lawDict)
        );

        if (laws_.size() == 1)
        {
            lawCells_[lawI] = Foam::identity(mesh_.nCells());
        }
        else // more than one material law
        {
            // Look up cell zone of the same name as the law
            const label zoneID = mesh_.cellZones().findZoneID(lawName);

            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "CellZone " << lawName << " not found"
                    << "When defining more than one mechanical constitutive "
                    << "law, each cell must belong to exactly one cellZone "
                    << "with the same name as the law."
                    << exit(FatalError);
            }

            lawCells_[lawI] = mesh_.cellZones()[zoneID];
        }

        // Check that each cell appears in only one cell zone
        forAll(lawCells_[lawI], i)
        {
            const label cellI = lawCells_[lawI][i];

            if (cellToLaw[cellI] != -1)
            {
                FatalErrorInFunction
                    << "Cell " << cellI
                    << " appears in more than one mechanical constitutive law: "
                    << lawNames[cellToLaw[cellI]] << " and " << lawName
                    << exit(FatalError);
            }

            cellToLaw[cellI] = lawI;
        }

        // State sized to this law's cells
        states_.set
        (
            lawI,
            new mechanicalConstitutiveLawState(lawCells_[lawI].size())
        );
    }

    // Check for any cells without a material law
    forAll(cellToLaw, cellI)
    {
        if (cellToLaw[cellI] == -1)
        {
            FatalErrorInFunction
                << "Cell " << cellI
                << " is not assigned to any mechanical constitutive law. "
                << "When defining more than one material, every cell must "
                << "belong to exactly one cellZone."
                << exit(FatalError);
        }
    }

    // Set lawBoundaryFaces
    forAll(lawNames, lawI)
    {
        lawBoundaryFaces_[lawI].resize(mesh.boundary().size());

        forAll(lawBoundaryFaces_, patchI)
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            labelHashSet curFaceSet;

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];

                if (cellToLaw[cellID] == lawI)
                {
                    curFaceSet.insert(faceI);
                }
            }

            lawBoundaryFaces_[lawI][patchI] = curFaceSet.toc();
        }

        // Initialise boundary patch states
        boundaryStates_[lawI].setSize(mesh_.boundary().size());
        forAll(mesh_.boundary(), patchI)
        {
            const label nFaces = lawBoundaryFaces_[lawI][patchI].size();

            boundaryStates_[lawI].set
            (
                patchI,
                new mechanicalConstitutiveLawState(nFaces)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::~mechanicalConstitutiveLawManager()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const volTensorField& gradD,
    const volTensorField& gradD0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    forAll(laws_, lawI)
    {
        const labelList& cells = lawCells_[lawI];

        // Update the internal field
        {
            // "View" into the gradD for this material => does not copy data
            const UIndirectList<tensor> gradDView
            (
                gradD.internalField(), cells
            );

            // "View" into the gradD0 for this material => does not copy data
            const UIndirectList<tensor> gradD0View
            (
                gradD0.internalField(), cells
            );

            // "View" into the stress for this material => does not copy data
            UIndirectList<symmTensor> stressView
            (
                stress.internalField(), cells
            );

            // Create wrapper for kinematic data: input to material law
            // This does not copy data
            smallStrainMechanicalConstitutiveLawKinematics kin
            (
                gradDView, gradD0View, dt
            );

            // Create wrapper for output
            if (scalarTangentPtr && tangentReq == tangentRequest::scalar)
            {
                // "View" into the tangent for this material => does not copy
                // data
                UIndirectList<scalar> tangentView
                (
                    scalarTangentPtr->internalField(),
                    cells
                );

                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView,
                    tangentView,
                    tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, states_[lawI], response);
            }
            else
            {
                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView,
                    tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, states_[lawI], response);
            }
        }

        // Update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        forAll(gradD.boundaryField(), patchI)
        {
            if (!gradD.boundaryField()[patchI].coupled())
            {
                // Select all faces on the patch for which the adjacent cell is
                // in this material
                const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                // "View" into the gradD for this material => does not copy data
                const UIndirectList<tensor> gradDView
                (
                    gradD.boundaryField()[patchI], faces
                );

                // "View" into the gradD0 for this material => does not copy data
                const UIndirectList<tensor> gradD0View
                (
                    gradD0.boundaryField()[patchI], faces
                );

                // "View" into the stress for this material => does not copy data
                UIndirectList<symmTensor> stressView
                (
                    stress.boundaryFieldRef()[patchI], faces
                );

                // Create wrapper for kinematic data: input to material law
                // This does not copy data
                smallStrainMechanicalConstitutiveLawKinematics kin
                (
                    gradDView, gradD0View, dt
                );

                // Create wrapper for output
                if (scalarTangentPtr && tangentReq == tangentRequest::scalar)
                {
                    // "View" into the tangent for this material => does not copy
                    // data
                    UIndirectList<scalar> tangentView
                    (
                        scalarTangentPtr->boundaryField()[patchI],
                        faces
                    );

                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView,
                        tangentView,
                        tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, boundaryStates_[lawI][patchI], response
                    );
                }
                else
                {
                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, boundaryStates_[lawI][patchI], response
                    );
                }
            }
        }
    }

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
(
    const volTensorField& F,
    const volTensorField& F0,
    const volScalarField& J,
    const volScalarField& J0,
    const volTensorField& Finv,
    const volTensorField& Finv0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    forAll(laws_, lawI)
    {
        const labelList& cells = lawCells_[lawI];

        // Update the internal field
        {
            // "View" into the F for this material => does not copy data
            const UIndirectList<tensor> FView(F.internalField(), cells);

            // "View" into the F0 for this material => does not copy data
            const UIndirectList<tensor> F0View(F0.internalField(), cells);

            // "View" into the J for this material => does not copy data
            const UIndirectList<scalar> JView(J.internalField(), cells);

            // "View" into the J0 for this material => does not copy data
            const UIndirectList<scalar> J0View(J0.internalField(), cells);

            // "View" into the Finv for this material => does not copy data
            const UIndirectList<tensor> FinvView(Finv.internalField(), cells);

            // "View" into the Finv0 for this material => does not copy data
            const UIndirectList<tensor> Finv0View(Finv0.internalField(), cells);

            // "View" into the stress for this material => does not copy data
            UIndirectList<symmTensor> stressView(stress.internalField(), cells);

            // Create wrapper for kinematic data: input to material law
            // This does not copy data
            finiteStrainMechanicalConstitutiveLawKinematics kin
            (
                FView, F0View, JView, J0View, FinvView, Finv0View, dt
            );

            // Create wrapper for output
            if (scalarTangentPtr && tangentReq == tangentRequest::scalar)
            {
                // "View" into the tangent for this material => does not copy
                // data
                UIndirectList<scalar> tangentView
                (
                    scalarTangentPtr->internalField(),
                    cells
                );

                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView,
                    tangentView,
                    tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, states_[lawI], response);
            }
            else
            {
                // Create wrapper for material: output from material law
                // This does not copy data
                mechanicalConstitutiveLawResponse response
                (
                    stressView,
                    tangentReq
                );

                // Update the material response, e.g. update the stress
                laws_[lawI].evaluate(kin, states_[lawI], response);
            }
        }

        // Update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        forAll(F.boundaryField(), patchI)
        {
            if (!F.boundaryField()[patchI].coupled())
            {
                // Select all faces on the patch for which the adjacent cell is
                // in this material
                const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                // "View" into the J for this material => does not copy data
                const UIndirectList<scalar> JView
                (
                    J.boundaryField()[patchI], faces
                );

                // "View" into the J0 for this material => does not copy data
                const UIndirectList<scalar> J0View
                (
                    J0.boundaryField()[patchI], faces
                );

                // "View" into the F for this material => does not copy data
                const UIndirectList<tensor> FView
                (
                    F.boundaryField()[patchI], faces
                );

                // "View" into the F0 for this material => does not copy data
                const UIndirectList<tensor> F0View
                (
                    F0.boundaryField()[patchI], faces
                );

                // "View" into the Finv for this material => does not copy data
                const UIndirectList<tensor> FinvView
                (
                    Finv.boundaryField()[patchI], faces
                );

                // "View" into the Finv0 for this material => does not copy data
                const UIndirectList<tensor> Finv0View
                (
                    Finv0.boundaryField()[patchI], faces
                );

                // "View" into the stress for this material => does not copy data
                UIndirectList<symmTensor> stressView
                (
                    stress.boundaryFieldRef()[patchI], faces
                );

                // Create wrapper for kinematic data: input to material law
                // This does not copy data
                finiteStrainMechanicalConstitutiveLawKinematics kin
                (
                    FView, F0View, JView, J0View, FinvView, Finv0View, dt
                );

                // Create wrapper for output
                if (scalarTangentPtr && tangentReq == tangentRequest::scalar)
                {
                    // "View" into the tangent for this material => does not copy
                    // data
                    UIndirectList<scalar> tangentView
                    (
                        scalarTangentPtr->boundaryField()[patchI],
                        faces
                    );

                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView,
                        tangentView,
                        tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, boundaryStates_[lawI][patchI], response
                    );
                }
                else
                {
                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, boundaryStates_[lawI][patchI], response
                    );
                }
            }
        }
    }

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();
}


// ************************************************************************* //
