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

#include "stabLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stabLeastSquaresVectors, 0);
}


const Foam::debug::tolerancesSwitch
Foam::stabLeastSquaresVectors::smallDotProdTol_
(
    "stabLeastSquaresSmallDotProdTol",
    0.1
);


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::stabLeastSquaresVectors::stabLeastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, stabLeastSquaresVectors>(mesh),
    pVectorsPtr_(nullptr),
    nVectorsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::stabLeastSquaresVectors::~stabLeastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stabLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "stabLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
    const volVectorField& C = mesh().C();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh().nCells(), symmTensor::zero);

    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const vector d = C[nei] - C[own];
        const symmTensor wdd = (1.0/magSqr(d))*sqr(d);

        dd[own] += wdd;
        dd[nei] += wdd;
    }

    forAll (lsP.boundaryField(), patchI)
    {
        const fvPatch& p = mesh().boundary()[patchI];
        const unallocLabelList& fc = p.faceCells();

        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        const vectorField pd = p.delta();

        forAll (pd, pFaceI)
        {
            const vector& d = pd[pFaceI];
            dd[fc[pFaceI]] += (1.0/magSqr(d))*sqr(d);
        }
    }

    // For easy access of neighbour coupled patch field needed for
    // lsN vectors on implicitly coupled boundaries. VV, 18/June/2014
    volSymmTensorField volInvDd
    (
        IOobject
        (
            "volInvDd",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero),
        zeroGradientFvPatchScalarField::typeName
    );
    symmTensorField& invDd = volInvDd.internalField();


    // Invert least squares matrix using Householder transformations to avoid
    // badly posed cells
    invDd = inv(dd); // ZT: as in extendedLeastSquares
    // invDd = hinv(dd);

    // Evaluate coupled to exchange coupled neighbour field data
    // across coupled boundaries.  HJ, 18/Mar/2015
    volInvDd.boundaryField().evaluateCoupled();

    // Revisit all faces and calculate the lsP and lsN vectors
    vectorField& lsPIn = lsP.internalField();
    vectorField& lsNIn = lsN.internalField();

    // Least squares vectors on internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const vector d = C[nei] - C[own];
        const scalar magSfByMagSqrd = 1.0/magSqr(d);

        lsPIn[faceI] = magSfByMagSqrd*(invDd[own] & d);
        lsNIn[faceI] = -magSfByMagSqrd*(invDd[nei] & d);
    }

    // Get inverse dd boundary field to grab the neighbouring side eventually
    const volSymmTensorField::GeometricBoundaryField& volInvDdb =
        volInvDd.boundaryField();

    // Least squares vectors on boundary faces
    forAll (lsP.boundaryField(), patchI)
    {
        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchI];
        fvsPatchVectorField& patchLsN = lsN.boundaryField()[patchI];
        const fvPatch& p = mesh().boundary()[patchI];
        const unallocLabelList& fc = p.faceCells();

        const fvPatchSymmTensorField& volInvDdp = volInvDdb[patchI];

        // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
        const vectorField pd = p.delta();

        if (volInvDdp.coupled())
        {
            const symmTensorField invDdNei = volInvDdp.patchNeighbourField();

            forAll (pd, pFaceI)
            {
                const vector& d = pd[pFaceI];

                patchLsP[pFaceI] = (1.0/magSqr(d))*(invDd[fc[pFaceI]] & d);

                patchLsN[pFaceI] = - (1.0/magSqr(d))*(invDdNei[pFaceI] & d);
            }
        }
        else
        {
            forAll (pd, pFaceI)
            {
                const vector& d = pd[pFaceI];

                patchLsP[pFaceI] = (1.0/magSqr(d))*(invDd[fc[pFaceI]] & d);
            }
        }
    }

    // Coefficient sign correction on least squares vectors
    // The sign of the coefficient must be positive to ensure correct
    // behaviour in coupled systems.  This is checked by evaluating
    // dot-product of the P/N vectors with the face area vector.
    // If the dot-product is negative, cell is marked for use with the
    // Gauss gradient, which is unconditionally positive
    // HJ, 21/Apr/2015

    // First loop: detect cells with bad least squares vectors

    // Use Gauss Gradient field: set to 1 if Gauss is needed
    volScalarField useGaussGrad
    (
        IOobject
        (
            "useGaussGrad",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    const surfaceVectorField& Sf = mesh().Sf();
    const surfaceScalarField& w = mesh().weights();

    const vectorField& SfIn = Sf.internalField();
    const scalarField& wIn = w.internalField();
    scalarField& uggIn = useGaussGrad.internalField();

    // Check internal faces
    forAll (owner, faceI)
    {
        if
        (
            (lsPIn[faceI] & SfIn[faceI])/(mag(lsPIn[faceI])*mag(SfIn[faceI]))
          < smallDotProdTol_()
        )
        {
            // Least square points in the wrong direction for owner
            // Use Gauss gradient
            uggIn[owner[faceI]] = 1;
        }

        if
        (
            (lsNIn[faceI] & SfIn[faceI])/(mag(lsNIn[faceI])*mag(SfIn[faceI]))
          > -smallDotProdTol_()
        )
        {
            // Least square points in the wrong direction for owner.
            // Note that Sf points into neighbour cell
            // Use Gauss gradient
            uggIn[neighbour[faceI]] = 1;
        }
    }

    forAll (lsP.boundaryField(), patchI)
    {
        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchI];
        const vectorField& pSf = Sf.boundaryField()[patchI];
        const fvPatch& p = mesh().boundary()[patchI];
        const unallocLabelList& fc = p.faceCells();

        // Same check for coupled and uncoupled
        forAll (patchLsP, pFaceI)
        {
            if
            (
                (patchLsP[pFaceI] & pSf[pFaceI])/
                (mag(patchLsP[pFaceI])*mag(pSf[pFaceI]))
              < smallDotProdTol_
            )
            {
                uggIn[fc[pFaceI]] = 1;
            }
        }
    }

    // Update boundary conditions for coupled boundaries.  This synchronises
    // the Gauss grad indication field
    useGaussGrad.boundaryField().evaluateCoupled();

    // Replace least square vectors with weighted Gauss gradient vectors
    // for marked cells

    // Prepare cell volumes with parallel communications
    volScalarField cellV
    (
        IOobject
        (
            "cellV",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    cellV.internalField() = mesh().V();

    // Evaluate coupled to exchange coupled neighbour field data
    // across coupled boundaries.  HJ, 18/Mar/2015
    cellV.boundaryField().evaluateCoupled();

    const scalarField& cellVIn = cellV.internalField();

    // Internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        if (uggIn[own] > SMALL)
        {
            // Gauss gradient for owner cell
            lsPIn[faceI] = (1 - wIn[faceI])*SfIn[faceI]/cellVIn[own];
        }

        if (uggIn[nei] > SMALL)
        {
            // Gauss gradient for neighbour cell
            lsNIn[faceI] = -wIn[faceI]*SfIn[faceI]/cellVIn[nei];
        }
    }

    // Get boundary field for cell volumes to fetch the other side eventually
    const volScalarField::GeometricBoundaryField& cellVb =
        cellV.boundaryField();

    // Boundary faces
    forAll (lsP.boundaryField(), patchI)
    {
        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchI];
        fvsPatchVectorField& patchLsN = lsN.boundaryField()[patchI];
        const fvPatch& p = mesh().boundary()[patchI];
        const unallocLabelList& fc = p.faceCells();

        const fvsPatchScalarField& pw = w.boundaryField()[patchI];
        const fvPatchScalarField& cellVp = cellVb[patchI];
        const vectorField& pSf = Sf.boundaryField()[patchI];

        // Get indicator for local side
        const fvPatchScalarField& ugg = useGaussGrad.boundaryField()[patchI];
        const scalarField pUgg = ugg.patchInternalField();

        if (cellVp.coupled())
        {
            const scalarField cellVInNei = cellVp.patchNeighbourField();

            // Get indicator for neighbour side
            const scalarField nUgg = ugg.patchNeighbourField();

            forAll (pUgg, pFaceI)
            {
                if (pUgg[pFaceI] > SMALL)
                {
                    // Gauss grad for owner cell
                    patchLsP[pFaceI] =
                        (1 - pw[pFaceI])*pSf[pFaceI]/cellVIn[fc[pFaceI]];
                }

                if (nUgg[pFaceI] > SMALL)
                {
                    // Gauss gradient for neighbour cell
                    patchLsN[pFaceI] =
                        -pw[pFaceI]*pSf[pFaceI]/cellVInNei[pFaceI];
                }
            }
        }
        else
        {
            forAll (pUgg, pFaceI)
            {
                if (pUgg[pFaceI] > SMALL)
                {
                    // Gauss grad for owner cell
                    patchLsP[pFaceI] =
                        (1 - pw[pFaceI])*pSf[pFaceI]/cellVIn[fc[pFaceI]];
                }
            }
        }
    }

    if (debug)
    {
        Info<< "stabLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField& Foam::stabLeastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::stabLeastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::stabLeastSquaresVectors::movePoints() const
{
    if (debug)
    {
        InfoIn("bool stabLeastSquaresVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}


bool Foam::stabLeastSquaresVectors::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn
        (
            "bool stabLeastSquaresVectors::updateMesh(const mapPolyMesh&) const"
        )   << "Clearing least square data" << endl;
    }

    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}

// ************************************************************************* //
