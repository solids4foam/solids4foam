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

\*---------------------------------------------------------------------------*/

#include "skewCorrectedLinearDC.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::skewCorrectedLinearDC<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
        gradVf = gradScheme_().grad(vf);


    // Since grad is used on coupled boundaries, correctBoundaryConditions
    // needs to be called.  HJ, 1/Nov/2012
    gradVf.correctBoundaryConditions();


    // Correction vectors

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

    const vectorField& Sf = mesh.Sf().internalField();
    const scalarField& magSf = mesh.magSf().internalField();

    vectorField nf = Sf/magSf;

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf&kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf&kNI)/sqr(magSf);

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
        }
    }


    // Perform skew-corrected interpolation

    Field<Type>& sfI = sf.internalField();
    const Field<Type>& vfI = vf.internalField();

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
    const scalarField& wI = w.internalField();

    forAll(sfI, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        Type vfOwn = vfI[own] + (kPI[faceI] & gradVf[own]);
        Type vfNei = vfI[nei] + (kNI[faceI] & gradVf[nei]);

        sfI[faceI] = wI[faceI]*(vfOwn - vfNei) + vfNei;
    }

    forAll (vf.boundaryField(), patchI)
    {
        const scalarField& pw =
            w.boundaryField()[patchI];

        if (vf.boundaryField()[patchI].coupled())
        {
            Field<Type> vfOwn =
                vf.boundaryField()[patchI].patchInternalField()
              + (
                    kP.boundaryField()[patchI]
                  & gradVf.boundaryField()[patchI].patchInternalField()
                );

            Field<Type> vfNei =
                vf.boundaryField()[patchI].patchNeighbourField()
              + (
                    kN.boundaryField()[patchI]
                  & gradVf.boundaryField()[patchI].patchNeighbourField()
                );

            sf.boundaryField()[patchI] = pw*(vfOwn - vfNei) + vfNei;
        }
        else
        {
            vf.boundaryField()[patchI].patchInterpolate(sf, pw);
        }
    }

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::skewCorrectedLinearDC<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Note: Correction is calculated by assembling the complete interpolation
    // including extrapolated gradient contribution and subtracting the
    // implicit contribution.  HJ, 27/Mar/2010
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "skewCorrectedLinearDCCorrection(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            skewCorrectedLinearDC<Type>::interpolate(vf)
          - upwind<Type>::interpolate
//           - surfaceInterpolationScheme<Type>::interpolate
            (
                vf,
                this->weights()
            )
        )
    );

    return tsfCorr;
}


namespace Foam
{
    //makelimitedSurfaceInterpolationScheme(skewCorrectedLinearDC)
    makelimitedSurfaceInterpolationTypeScheme(skewCorrectedLinearDC, scalar)
    makelimitedSurfaceInterpolationTypeScheme(skewCorrectedLinearDC, vector)
}

// ************************************************************************* //
