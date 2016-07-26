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

Description
    Zoned multi-material mechanical controlled by a material indicator field.

\*---------------------------------------------------------------------------*/

#include "multiMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterial, 0);
    addToRunTimeSelectionTable(mechanicalLaw, multiMaterial, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multiMaterial::indicator
(
    const label i
) const
{
    const scalarField& mat = materials_.internalField();

    tmp<scalarField> tresult(new scalarField(mat.size(), 0.0));
    scalarField& result = tresult();

    forAll(mat, matI)
    {
        if (mat[matI] > i - SMALL && mat[matI] < i + 1 - SMALL)
        {
            result[matI] = 1.0;
        }
    }

    return tresult;
}


Foam::scalar
Foam::multiMaterial::indicator(const label index, const label cellID) const
{
    const scalar mat = materials_.internalField()[cellID];
    scalar result = 0.0;

    if (mat > index - SMALL && mat < index + 1 - SMALL)
      {
    result = 1.0;
      }

    return result;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiMaterial::multiMaterial
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    PtrList<mechanicalLaw>(),
    materials_
    (
        IOobject
        (
            "materials",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    PtrList<mechanicalLaw>& laws = *this;

    PtrList<entry> lawEntries(dict.lookup("laws"));
    laws.setSize(lawEntries.size());

    forAll(laws, lawI)
    {
        laws.set
        (
            lawI,
            mechanicalLaw::New
            (
                lawEntries[lawI].keyword(),
                mesh,
                lawEntries[lawI].dict()
            )
        );

        laws[lawI].setMaterialIndex(lawI);
    }

    if
    (
        min(materials_).value() < 0
     || max(materials_).value() > (laws.size() - SMALL)
    )
    {
        FatalErrorIn
        (
            "multiMaterial::multiMaterial\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Invalid definition of material indicator field.  "
            << "There are " << laws.size() << " defined in " << dict.name()
            << nl << "but the max/min material index lies outside this range"
            << nl << "Max material index is " << min(materials_).value() << nl
            << "Min material index is " << max(materials_).value() << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiMaterial::~multiMaterial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiMaterial::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroRho", dimDensity, 0),
            calculatedFvPatchScalarField::typeName
            //zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].rho()().internalField();
    }

    forAll(result.boundaryField(), patchI)
    {
        if
        (
            !result.boundaryField()[patchI].coupled()
         && result.boundaryField()[patchI].type() != "empty"
        )
        {
            result.boundaryField()[patchI] =
                result.boundaryField()[patchI].patchInternalField();
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterial::impK() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimForce/dimArea, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].impK()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterial::E() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroE", dimForce/dimArea, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].E()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField>
Foam::multiMaterial::E(const volScalarField& epsEq) const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroE", dimForce/dimArea, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].E(epsEq)().internalField();
    }

//     forAll(result.boundaryField(),patchI)
//     {
//         forAll(laws, lawI)
//         {
//             result.boundaryField()[patchI] +=
//                 indicator(lawI)().boundaryField()[patchI]
//                 *laws[lawI].E(t)().boundaryField()[patchI];
//         }
//     }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterial::nu() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroE", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].nu()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volDiagTensorField> Foam::multiMaterial::K() const
{
    tmp<volDiagTensorField> tresult
    (
        new volDiagTensorField
        (
            IOobject
            (
                "K",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedDiagTensor("zeroK", dimForce/dimArea, diagTensor::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volDiagTensorField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
      {
    result.internalField() +=
      indicator(lawI)*laws[lawI].K()().internalField();
      }

    result.correctBoundaryConditions();

    return tresult;
}

Foam::tmp<Foam::volSymmTensor4thOrderField> Foam::multiMaterial::C() const
{
    tmp<volSymmTensor4thOrderField> tresult
    (
        new volSymmTensor4thOrderField
        (
            IOobject
            (
                "C",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor4thOrder
            ("zeroC", dimForce/dimArea, symmTensor4thOrder::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volSymmTensor4thOrderField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result.internalField() +=
      indicator(lawI)*laws[lawI].C()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


void Foam::multiMaterial::correct()
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].correct();
    }
}


void Foam::multiMaterial::correct(volSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].correct(sigma);
    }
}


void Foam::multiMaterial::correct(volSymmTensorField& sigma, const int flag)
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].correct(sigma, flag);
    }
}


void Foam::multiMaterial::correct(surfaceSymmTensorField& sigma)
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].correct(sigma);
    }
}


Foam::scalar Foam::multiMaterial::residual()
{
    PtrList<mechanicalLaw>& laws = *this;

    scalar maxRes = 0.0;

    forAll(laws, lawI)
    {
        maxRes = max(maxRes, laws[lawI].residual());
    }

    return maxRes;
}


void Foam::multiMaterial::updateYieldStress()
{
    PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        laws[lawI].updateYieldStress();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::multiMaterial::plasticDissipationRate() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "PlasticDissipationRate",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroE", dimForce/(dimArea*dimTime), 0)
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<mechanicalLaw>& laws = *this;

    forAll(laws, lawI)
    {
        result += laws[lawI].curMaterial()*laws[lawI].plasticDissipationRate();
    }

    result.correctBoundaryConditions();

    return tresult;
}


// ************************************************************************* //
