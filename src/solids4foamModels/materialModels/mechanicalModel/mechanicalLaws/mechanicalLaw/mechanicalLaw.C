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

Foam::tmp<Foam::volScalarField> Foam::mechanicalLaw::impK() const
{
    const volScalarField E = this->E();
    const volScalarField nu = this->nu();

    const volScalarField mu = E/(2.0*(1.0 + nu));

    if
    (
        mesh_.lookupObject<IOdictionary>
        (
            "mechanicalProperties"
        ).lookup("planeStress")
    )
    {
        const volScalarField lambda = nu*E/((1.0 + nu)*(1.0 - nu));

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "impK",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                2.0*mu + lambda
            )
        );
    }
    else
    {
        const volScalarField lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

            return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "impK",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                2.0*mu + lambda
            )
        );
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalLaw::impKf() const
{
    return fvc::interpolate(impK());
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalLaw::gammaInf() const
{
    notImplemented(type() + "::gammaInf()");

    return E();
}


Foam::PtrList<Foam::volScalarField> Foam::mechanicalLaw::gamma() const
{
    notImplemented(type() + "::gamma()");

    return PtrList<volScalarField>(0);
}


Foam::PtrList<Foam::volScalarField> Foam::mechanicalLaw::tau() const
{
    notImplemented(type() + "::tau()");

    return PtrList<volScalarField>(0);
}

Foam::tmp<Foam::volDiagTensorField> Foam::mechanicalLaw::K() const
{
  volScalarField J =
      ( (1 - nu()*nu() - nu()*nu() - nu()*nu() - 2*nu()*nu()*nu())
        /(E()*E()*E()) );
  volScalarField A11 = ( (1 - nu()*nu())/(J*E()*E()) );
  volScalarField A22 = ( (1 - nu()*nu())/(J*E()*E()) );
  volScalarField A33 = ( (1 - nu()*nu())/(J*E()*E()) );

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
        dimensionedDiagTensor("K", A11.dimensions(), diagTensor::zero)
        )
    );

  volDiagTensorField& result = tresult();

  result.replace(diagTensor::XX, A11);
  result.replace(diagTensor::YY, A22);
  result.replace(diagTensor::ZZ, A33);

  tresult().correctBoundaryConditions();

  return tresult;
}


Foam::tmp<Foam::volSymmTensor4thOrderField> Foam::mechanicalLaw::C() const
{
  volScalarField twoMu = 2.0*E()/(2.0*(1+nu()));
  volScalarField lambda = nu()*E()/((1.0 + nu())*(1.0 - 2.0*nu()));
  volScalarField twoMuLambda = twoMu + lambda;

  // symmTensor4thOrder C (
  //            2*mu + lambda, lambda, lambda,
  //            2*mu + lambda, lambda,
  //            2*mu + lambda,
  //            2*mu,
  //            2*mu,
  //            2*mu
  //            );

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
            ("C", twoMu.dimensions(), symmTensor4thOrder::zero)
        )
    );

  volSymmTensor4thOrderField& result = tresult();

  result.replace(symmTensor4thOrder::XXXX, twoMuLambda);
  result.replace(symmTensor4thOrder::XXYY, lambda);
  result.replace(symmTensor4thOrder::XXZZ, lambda);

  result.replace(symmTensor4thOrder::YYYY, twoMuLambda);
  result.replace(symmTensor4thOrder::YYZZ, lambda);

  result.replace(symmTensor4thOrder::ZZZZ, twoMuLambda);

  result.replace(symmTensor4thOrder::XYXY, twoMu);
  result.replace(symmTensor4thOrder::YZYZ, twoMu);
  result.replace(symmTensor4thOrder::ZXZX, twoMu);

  tresult().correctBoundaryConditions();

  return tresult;
}


tmp<volScalarField> Foam::mechanicalLaw::bulkModulus() const
{
    notImplemented(type() + "::bulkModulus()");

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
              "bulkModulus",
              mesh().time().timeName(),
              mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    return tresult;
}


tmp<volScalarField> Foam::mechanicalLaw::viscosity() const
{
    notImplemented(type() + "::viscosity()");

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
              "viscosity",
              mesh().time().timeName(),
              mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    return tresult;
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
