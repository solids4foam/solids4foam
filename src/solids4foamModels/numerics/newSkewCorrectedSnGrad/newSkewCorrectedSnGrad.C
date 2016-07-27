/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Simple central-difference snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "newSkewCorrectedSnGrad.H"
#include "fv.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "correctedSnGrad.H"
#include "uncorrectedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
newSkewCorrectedSnGrad<Type>::~newSkewCorrectedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
newSkewCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tssf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf();

    ssf = dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero);

//     typedef typename pTraits<Type>::cmptType cType;

    typedef typename
        outerProduct<vector, typename pTraits<Type>::cmptType>::type 
        CmptGradType;

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const vectorField& Sf = mesh.Sf().internalField();
    const scalarField& magSf = mesh.magSf().internalField();

    vectorField nf = Sf/magSf;

    const vectorField& Cf = mesh.Cf().internalField();
    const vectorField& C = mesh.C().internalField();

    const scalarField& deltaCoeffs = 
        mesh.deltaCoeffs().internalField();

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
    vectorField& kPI = kP.internalField();

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
    vectorField& kNI = kN.internalField();

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf&kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf&kNI)/sqr(magSf);

//     vectorField delta = 
//         Cf 
//       - (vectorField(C, neighbour) + kN + vectorField(C, owner) + kP)/2.0;

//     kPI += delta;
//     kNI += delta;

    forAll(kP.boundaryField(), patchI)
    {
        if (kP.boundaryField()[patchI].coupled())
        {
            kP.boundaryField()[patchI] = 
                mesh.boundary()[patchI].Cf()
              - mesh.boundary()[patchI].Cn();

            kP.boundaryField()[patchI] -= 
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kP.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

            kN.boundaryField()[patchI] = 
                mesh.Cf().boundaryField()[patchI]
              - (
                    mesh.boundary()[patchI].Cn()
                  + mesh.boundary()[patchI].delta()
                );

            kN.boundaryField()[patchI] -= 
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kN.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

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

//             kP.boundaryField()[patchI] += delta;
//             kN.boundaryField()[patchI] += delta;
        }
    }

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        GeometricField<CmptGradType, fvPatchField, volMesh> cmptGrad =
            gradScheme<typename pTraits<Type>::cmptType>::New
            (
                mesh,
                mesh.schemesDict().gradScheme(ssf.name())
            )()
           .grad(vf.component(cmpt));

        const Field<CmptGradType>& cmptGradI = cmptGrad.internalField();

        // Skewness and non-rothogonal correction
        {
            ssf.internalField().replace
            (
                cmpt,
                (
                    (kNI&Field<CmptGradType>(cmptGradI, neighbour))
                  - (kPI&Field<CmptGradType>(cmptGradI, owner))
                )
               *deltaCoeffs
            );
        }

        forAll(ssf.boundaryField(), patchI)
        {
            if (ssf.boundaryField()[patchI].coupled())
            {
                ssf.boundaryField()[patchI].replace
                (
                    cmpt,
                    (
                        (
                            kN.boundaryField()[patchI]
                          & cmptGrad.boundaryField()[patchI]
                           .patchNeighbourField()
                        )
                      - (
                            kP.boundaryField()[patchI]
                          & cmptGrad.boundaryField()[patchI]
                           .patchInternalField()
                        )
                    )
                   *mesh.deltaCoeffs().boundaryField()[patchI]
                );                
            }
        }
    }

    surfaceScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag
            (
                uncorrectedSnGrad<Type>::snGrad
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
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
    }

    ssf *= limiter;

    return tssf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
