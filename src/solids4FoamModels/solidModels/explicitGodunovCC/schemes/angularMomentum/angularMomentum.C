/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "angularMomentum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(angularMomentum, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

angularMomentum::angularMomentum
(
    const fvMesh& vm,
    const dictionary& dict
)
:
    mesh_(vm),
    rho_(dict.lookup("rho"))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
angularMomentum::~angularMomentum()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void angularMomentum::AMconservation
(
    GeometricField<vector, fvPatchField, volMesh>& rhsLm,
    GeometricField<vector, fvPatchField, volMesh>& rhsLm1,
    const GeometricField<vector, fvPatchField, volMesh>& rhsAm,
    const scalar& RKstage
) const
{
    const scalarField& V_ = mesh_.V();
    const objectRegistry& db = mesh_.thisDb();
    const volVectorField& x_ = db.lookupObject<volVectorField>("x");
    const volVectorField& lm_ = db.lookupObject<volVectorField>("lm");

    const dimensionedScalar deltaT
    (
        "deltaT",
        dimensionSet(0,0,1,0,0,0,0),
        db.time().deltaTValue()
    );

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf_x
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "x",
                x_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("x", x_.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> xAM = tvf_x();

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf_lm
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "lm",
                lm_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("lm", lm_.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> lmAM = tvf_lm();

    if (RKstage == 0)
    {
        xAM = x_.oldTime();
    }

    else if (RKstage == 1)
    {
        xAM = x_.oldTime() + (deltaT/2.0)*(lm_.oldTime()/rho_);
        lmAM = lm_.oldTime() + (deltaT*rhsLm1);
        xAM = xAM + ((deltaT*(lmAM/rho_))/2.0);
    }

    tensor K_LL = tensor::zero;
    tensor K_LB = tensor::zero;
    scalar K_BB = 0.0;
    vector R_L = vector::zero;

    forAll(mesh_.cells(), cellID)
    {
        K_LL +=
            V_[cellID]
           *((xAM[cellID] & xAM[cellID])*tensor::I - (xAM[cellID]*xAM[cellID]));

        K_LB += V_[cellID]
           *tensor(0, -xAM[cellID].z(), xAM[cellID].y(), xAM[cellID].z(), 0,
           -xAM[cellID].x(), -xAM[cellID].y(), xAM[cellID].x(), 0);

        K_BB += -V_[cellID];

        R_L +=
            (V_[cellID]*rhsAm[cellID])
          + ((V_[cellID]*rhsLm[cellID]) ^ xAM[cellID]);
    }

    if (Pstream::parRun())
    {
        reduce(K_LL, sumOp<tensor>());
        reduce(K_LB, sumOp<tensor>());
        reduce(K_BB, sumOp<scalar>());
        reduce(R_L, sumOp<vector>());
    }

    tensor LHS = K_LL - ((K_LB & K_LB)/K_BB);
    vector RHS = R_L;

    vector lambda = inv(LHS) & RHS;
    vector beta = (-K_LB & lambda)/K_BB;

    forAll(mesh_.cells(), cellID)
    {
        rhsLm[cellID] = rhsLm[cellID] + (lambda ^ xAM[cellID]) + beta;
    }

    if (RKstage == 0)
    {
        rhsLm1 = rhsLm;
    }

    tvf_x.clear();
    tvf_lm.clear();

    if (Pstream::parRun())
    {
        rhsLm.correctBoundaryConditions();
        rhsLm1.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void angularMomentum::printGlobalMomentum
(
    const GeometricField<vector, fvPatchField, volMesh>& lm,
    const GeometricField<vector, fvPatchField, volMesh>& x
) const
{

    vector lmG = vector::zero;
    vector amG = vector::zero;
    scalar vol = gSum(mesh_.V());

    forAll(mesh_.cells(), cellID)
    {
        lmG += lm[cellID]*mesh_.V()[cellID];
        amG += mesh_.V()[cellID]*(x[cellID] ^ lm[cellID]);
    }

    if (Pstream::parRun())
    {
        reduce(lmG, sumOp<vector>());
        reduce(amG, sumOp<vector>());
    }

    Info<< "\nPrinting global momentums ..." << nl
        << "Global linear momentum = " << lmG/vol << nl
        << "Global angular momentum = " << amG/vol << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //