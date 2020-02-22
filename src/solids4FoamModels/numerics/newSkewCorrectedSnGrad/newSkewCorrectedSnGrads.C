/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Simple central-difference snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "newSkewCorrectedSnGrad.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef OPENFOAMESIORFOUNDATION
makeSnGradScheme(newSkewCorrectedSnGrad)
#endif

namespace Foam
{

namespace fv
{

#ifndef OPENFOAMESIORFOUNDATION
makeSnGradScheme(newSkewCorrectedSnGrad)
#endif

template<>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >
newSkewCorrectedSnGrad<scalar>::correction
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<scalar, fvsPatchField, surfaceMesh>
    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > tssf
    (
        new GeometricField<scalar, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "snGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
#ifdef OPENFOAMESIORFOUNDATION
    GeometricField<scalar, fvsPatchField, surfaceMesh>& ssf = tssf.ref();
#else
    GeometricField<scalar, fvsPatchField, surfaceMesh>& ssf = tssf();
#endif

    ssf = dimensioned<scalar>("0", ssf.dimensions(), 0);

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

#ifdef OPENFOAMESIORFOUNDATION
    const vectorField& Sf = mesh.Sf().primitiveField();
//     const scalarField& magSf = mesh.magSf().primitiveField();
#else
    const vectorField& Sf = mesh.Sf().internalField();
//     const scalarField& magSf = mesh.magSf().internalField();
#endif

    vectorField nf = Sf/mag(Sf);

#ifdef OPENFOAMESIORFOUNDATION
    const vectorField& Cf = mesh.Cf().primitiveField();
    const vectorField& C = mesh.C().primitiveField();

    const scalarField& deltaCoeffs =
        mesh.deltaCoeffs().primitiveField();
#else
    const vectorField& Cf = mesh.Cf().internalField();
    const vectorField& C = mesh.C().internalField();

    const scalarField& deltaCoeffs =
        mesh.deltaCoeffs().internalField();
#endif

    surfaceVectorField kP
    (
        IOobject
        (
            "kP",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    );
#ifdef OPENFOAMESIORFOUNDATION
    vectorField& kPI = kP.primitiveFieldRef();
#else
    vectorField& kPI = kP.internalField();
#endif

    surfaceVectorField kN
    (
        IOobject
        (
            "kN",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    );
#ifdef OPENFOAMESIORFOUNDATION
    vectorField& kNI = kN.primitiveFieldRef();
#else
    vectorField& kNI = kN.internalField();
#endif

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf & kPI)/magSqr(Sf);
//     kPI -= Sf*(Sf&kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf & kNI)/magSqr(Sf);
//     kNI -= Sf*(Sf&kNI)/sqr(magSf);

//     vectorField delta =
//         Cf
//       - (vectorField(C, neighbour) + kN + vectorField(C, owner) + kP)/2.0;

//     kPI += delta;
//     kNI += delta;

    forAll(kP.boundaryField(), patchI)
    {
        if (kP.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            kP.boundaryFieldRef()[patchI] =
#else
            kP.boundaryField()[patchI] =
#endif
                mesh.boundary()[patchI].Cf()
              - mesh.boundary()[patchI].Cn();

#ifdef OPENFOAMESIORFOUNDATION
            kP.boundaryFieldRef()[patchI] -=
#else
            kP.boundaryField()[patchI] -=
#endif
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kP.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

#ifdef OPENFOAMESIORFOUNDATION
            kN.boundaryFieldRef()[patchI] =
#else
            kN.boundaryField()[patchI] =
#endif
                mesh.Cf().boundaryField()[patchI]
              - (
                    mesh.boundary()[patchI].Cn()
                  + mesh.boundary()[patchI].delta()
                );

#ifdef OPENFOAMESIORFOUNDATION
            kN.boundaryFieldRef()[patchI] -=
#else
            kN.boundaryField()[patchI] -=
#endif
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kN.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

//             Info << mesh.boundary()[patchI].name() << endl;
//             Info << 1.0/mesh.boundary()[patchI].deltaCoeffs() << endl;

//             vectorField delta =
//                 mesh.boundary()[patchI].Cf()
//               - (
//                     (
//                         mesh.boundary()[patchI].Cn()
//                       + mesh.boundary()[patchI].delta()
//                     )
//                   + kN.boundaryField()[patchI]
//                   + mesh.boundary()[patchI].Cn()
//                   + kP.boundaryField()[patchI]
//                 )/2.0;

// #ifdef OPENFOAMESIORFOUNDATION
//             kP.boundaryFieldRef()[patchI] += delta;
//             kN.boundaryFieldRef()[patchI] += delta;
// #else
//             kP.boundaryField()[patchI] += delta;
//             kN.boundaryField()[patchI] += delta;
// #endif
        }
    }

    const volVectorField* gradVfPtr(NULL);

    bool calcGradient = false;

    if ( mesh.foundObject<volVectorField>("grad(" + vf.name() + ")") )
    {
        const volVectorField& gradVf_ =
            mesh.lookupObject<volVectorField>("grad(" + vf.name() + ")");

        gradVfPtr = &gradVf_;
    }
    else
    {
        gradVfPtr =
            new volVectorField
            (
                gradScheme<scalar>::New
                (
                    mesh,
#ifdef OPENFOAMESIORFOUNDATION
                    mesh.gradScheme(ssf.name())
#else
                    mesh.schemesDict().gradScheme(ssf.name())
#endif
                )()
               .grad(vf)
            );

        calcGradient = true;
    }

    const volVectorField& gradVf = *gradVfPtr;

#ifdef OPENFOAMESIORFOUNDATION
    ssf.primitiveFieldRef() =
#else
    ssf.internalField() =
#endif
        (
            (kNI & Field<vector>(gradVf, neighbour))
          - (kPI & Field<vector>(gradVf, owner))
        )
       *deltaCoeffs;

    forAll(ssf.boundaryField(), patchI)
    {
        if (ssf.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            ssf.boundaryFieldRef()[patchI] =
#else
            ssf.boundaryField()[patchI] =
#endif
            (
                (
                    (
                        kN.boundaryField()[patchI]
                      & gradVf.boundaryField()[patchI].patchNeighbourField()
                    )
                  - (
                        kP.boundaryField()[patchI]
                      & gradVf.boundaryField()[patchI].patchInternalField()
                    )
                )
               *mesh.deltaCoeffs().boundaryField()[patchI]
            );
        }
    }

    surfaceScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag
            (
                uncorrectedSnGrad<scalar>::snGrad
                (
                    vf,
                    this->deltaCoeffs(vf),
                    "orthSnGrad"
                )
              + ssf
            )
           /(
                (1 - limitCoeff_)*mag(ssf)
              + dimensionedScalar("small", ssf.dimensions(), SMALL)
            ),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    if (fv::debug)
    {
#ifdef OPENFOAMESIORFOUNDATION
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.primitiveField())
            << " max: "<< max(limiter.primitiveField())
            << " avg: " << average(limiter.primitiveField()) << endl;
#else
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
#endif
    }

    ssf *= limiter;

    if (calcGradient)
    {
        deleteDemandDrivenData(gradVfPtr);
    }

    return tssf;
}


template<>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> >
newSkewCorrectedSnGrad<vector>::correction
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<vector, fvsPatchField, surfaceMesh>
    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tssf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "snGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
#ifdef OPENFOAMESIORFOUNDATION
    GeometricField<vector, fvsPatchField, surfaceMesh>& ssf = tssf.ref();
#else
    GeometricField<vector, fvsPatchField, surfaceMesh>& ssf = tssf();
#endif

    ssf = dimensioned<vector>("0", ssf.dimensions(), vector::zero);

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

#ifdef OPENFOAMESIORFOUNDATION
    const vectorField& Sf = mesh.Sf().primitiveField();
//     const scalarField& magSf = mesh.magSf().primitiveField();
#else
    const vectorField& Sf = mesh.Sf().internalField();
//     const scalarField& magSf = mesh.magSf().internalField();
#endif

    vectorField nf = Sf/mag(Sf);

#ifdef OPENFOAMESIORFOUNDATION
    const vectorField& Cf = mesh.Cf().primitiveField();
    const vectorField& C = mesh.C().primitiveField();

    const scalarField& deltaCoeffs =
        mesh.deltaCoeffs().primitiveField();
#else
    const vectorField& Cf = mesh.Cf().internalField();
    const vectorField& C = mesh.C().internalField();

    const scalarField& deltaCoeffs =
        mesh.deltaCoeffs().internalField();
#endif

    surfaceVectorField kP
    (
        IOobject
        (
            "kP",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    );
#ifdef OPENFOAMESIORFOUNDATION
    vectorField& kPI = kP.primitiveFieldRef();
#else
    vectorField& kPI = kP.internalField();
#endif

    surfaceVectorField kN
    (
        IOobject
        (
            "kN",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    );
#ifdef OPENFOAMESIORFOUNDATION
    vectorField& kNI = kN.primitiveFieldRef();
#else
    vectorField& kNI = kN.internalField();
#endif

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf & kPI)/magSqr(Sf);
//     kPI -= Sf*(Sf & kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf & kNI)/magSqr(Sf);
//     kNI -= Sf*(Sf & kNI)/sqr(magSf);

    // non-uniformity correction
    if (false)
    {
        vectorField delta =
            Cf
          - (vectorField(C, neighbour) + kN + vectorField(C, owner) + kP)/2.0;

        kPI += delta;
        kNI += delta;
    }

    forAll(kP.boundaryField(), patchI)
    {
        if (kP.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            kP.boundaryFieldRef()[patchI] =
#else
            kP.boundaryField()[patchI] =
#endif
                mesh.boundary()[patchI].Cf()
              - mesh.boundary()[patchI].Cn();

#ifdef OPENFOAMESIORFOUNDATION
            kP.boundaryFieldRef()[patchI] -=
#else
            kP.boundaryField()[patchI] -=
#endif
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kP.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

#ifdef OPENFOAMESIORFOUNDATION
            kN.boundaryFieldRef()[patchI] =
#else
            kN.boundaryField()[patchI] =
#endif
                mesh.Cf().boundaryField()[patchI]
              - (
                    mesh.boundary()[patchI].Cn()
                  + mesh.boundary()[patchI].delta()
                );

#ifdef OPENFOAMESIORFOUNDATION
            kN.boundaryFieldRef()[patchI] -=
#else
            kN.boundaryField()[patchI] -=
#endif
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kN.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

//             Info << mesh.boundary()[patchI].name() << endl;
//             Info << 1.0/mesh.boundary()[patchI].deltaCoeffs() << endl;

//             // Ununiformity corr
//             vectorField delta =
//                 mesh.boundary()[patchI].Cf()
//               - (
//                     (
//                         mesh.boundary()[patchI].Cn()
//                       + mesh.boundary()[patchI].delta()
//                     )
//                   + kN.boundaryField()[patchI]
//                   + mesh.boundary()[patchI].Cn()
//                   + kP.boundaryField()[patchI]
//                 )/2.0;

// #ifdef OPENFOAMESIORFOUNDATION
//             kP.boundaryFieldRef()[patchI] += delta;
//             kN.boundaryFieldRef()[patchI] += delta;
// #else
//             kP.boundaryField()[patchI] += delta;
//             kN.boundaryField()[patchI] += delta;
// #endif
        }
    }

    const volTensorField* gradVfPtr(NULL);

    bool calcGradient = false;

    if ( mesh.foundObject<volTensorField>("grad(" + vf.name() + ")") )
    {
        const volTensorField& gradVf_ =
            mesh.lookupObject<volTensorField>("grad(" + vf.name() + ")");

        gradVfPtr = &gradVf_;
    }
    else
    {
        gradVfPtr =
            new volTensorField
            (
                gradScheme<vector>::New
                (
                    mesh,
#ifdef OPENFOAMESIORFOUNDATION
                    mesh.gradScheme(ssf.name())
#else
                    mesh.schemesDict().gradScheme(ssf.name())
#endif
                )()
               .grad(vf)
            );

        calcGradient = true;
    }

    const volTensorField& gradVf = *gradVfPtr;

#ifdef OPENFOAMESIORFOUNDATION
    ssf.primitiveFieldRef() =
#else
    ssf.internalField() =
#endif
        (
            (kNI & Field<tensor>(gradVf, neighbour))
          - (kPI & Field<tensor>(gradVf, owner))
        )
       *deltaCoeffs;

    forAll(ssf.boundaryField(), patchI)
    {
        if (ssf.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            ssf.boundaryFieldRef()[patchI] =
#else
            ssf.boundaryField()[patchI] =
#endif
            (
                (
                    (
                        kN.boundaryField()[patchI]
                      & gradVf.boundaryField()[patchI].patchNeighbourField()
                    )
                  - (
                        kP.boundaryField()[patchI]
                      & gradVf.boundaryField()[patchI].patchInternalField()
                    )
                )
               *mesh.deltaCoeffs().boundaryField()[patchI]
            );
        }
    }

    surfaceScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag
            (
                uncorrectedSnGrad<vector>::snGrad
                (
                    vf,
                    this->deltaCoeffs(vf),
                    "orthSnGrad"
                )
              + ssf
            )
           /(
                (1 - limitCoeff_)*mag(ssf)
              + dimensionedScalar("small", ssf.dimensions(), SMALL)
            ),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    if (fv::debug)
    {
#ifdef OPENFOAMESIORFOUNDATION
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.primitiveField())
            << " max: "<< max(limiter.primitiveField())
            << " avg: " << average(limiter.primitiveField()) << endl;
#else
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
#endif
    }

    ssf *= limiter;

    if (calcGradient)
    {
        deleteDemandDrivenData(gradVfPtr);
    }

    return tssf;
}


}
}

// ************************************************************************* //
