/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    mechanicalLaw

Description
    Material mechanical for solids.

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanicalLaw, 0);
defineRunTimeSelectionTable(mechanicalLaw, dictionary);


// * * * * * * * * * * * Private Member functions * * * * * * * * * * * * * * //

void mechanicalLaw::calcCurMaterial() const
{
    if (curMaterialPtr_)
    {
        FatalErrorIn
        (
            "const volScalarField& plasticityStressReturn::curMaterial() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    if (mesh.objectRegistry::foundObject<volScalarField>("materials"))
    {
        curMaterialPtr_ =
            new volScalarField
            (
                IOobject
                (
                    "curMaterial",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
            );

        volScalarField& curMaterial = *curMaterialPtr_;

        scalarField& curMaterialI = curMaterial.internalField();

        // Check the current material index has been set
        if (curMatIndex_ < 0)
        {
            FatalErrorIn("void mechanicalLaw::calcCurMaterial() const")
                << "The current material index has not been set"
                << abort(FatalError);
        }

        // Lookup materials index field
        const volScalarField& materials =
            mesh.objectRegistry::lookupObject<volScalarField>("materials");
        const scalarField& materialsI = materials.internalField();

        forAll(materialsI, cellI)
        {
            if (label(materialsI[cellI]) == curMatIndex_)
            {
                curMaterialI[cellI] = 1.0;
            }
        }

        forAll(materials.boundaryField(), patchI)
        {
            if
            (
                !materials.boundaryField()[patchI].coupled()
             && mesh.boundaryMesh().type() != emptyPolyPatch::typeName
            )
            {
                scalarField& pCurMaterial =
                    curMaterial.boundaryField()[patchI];
                const labelList& faceCells =
                    mesh.boundaryMesh()[patchI].faceCells();

                forAll(pCurMaterial, faceI)
                {
                    if (label(materialsI[faceCells[faceI]]) == curMatIndex_)
                    {
                        pCurMaterial[faceI] = 1.0;
                    }
                }
            }
        }

        curMaterial.correctBoundaryConditions();
    }
    else
    {
        curMaterialPtr_ =
            new volScalarField
            (
                IOobject
                (
                    "curMaterial",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("1.0", dimless, 1.0),
                zeroGradientFvPatchScalarField::typeName
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mechanicalLaw::mechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    curMaterialPtr_(NULL),
    curMatIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalLaw::impKf() const
{
    return fvc::interpolate(impK());
}


void Foam::mechanicalLaw::correct(surfaceSymmTensorField&)
{
    FatalErrorIn(type() + "::correct(surfaceSymmTensorField&)")
        << "The correct(surfaceSymmTensorField&) function is not defined"
        << " for the current mechanical law" << abort(FatalError);
}


const Foam::volScalarField& Foam::mechanicalLaw::curMaterial() const
{
    if (!curMaterialPtr_)
    {
        calcCurMaterial();
    }

    return *curMaterialPtr_;
}


Foam::volScalarField& Foam::mechanicalLaw::curMaterial()
{
    if (!curMaterialPtr_)
    {
        calcCurMaterial();
    }

    return *curMaterialPtr_;
}


Foam::scalar Foam::mechanicalLaw::residual()
{
    // Default to zero; this can be overwritten by any derived mechanical law
    // For example in a elasto-plastic law, the plastic increment may need to be
    // checked
    return 0.0;
}


} // End namespace Foam

// ************************************************************************* //
