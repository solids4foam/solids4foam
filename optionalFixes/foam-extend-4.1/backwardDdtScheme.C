/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "backwardDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
// #include "slipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar backwardDdtScheme<Type>::deltaT_() const
{
    return mesh().time().deltaT().value();
}


template<class Type>
scalar backwardDdtScheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


template<class Type>
template<class GeoField>
scalar backwardDdtScheme<Type>::deltaT0_(const GeoField& vf) const
{
    // Bug fix, Zeljko Tukovic: solver with outer iterations over a time-step
    // HJ, 12/Feb/2010
    if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
backwardDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt().internalField() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().V0() - coefft00*mesh().V00())/mesh().V()
        );

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
backwardDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime().internalField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().internalField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
backwardDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.internalField() -
                    (
                        coefft0*vf.oldTime().internalField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().internalField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
backwardDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*rho.internalField()*vf.internalField() -
                    (
                        coefft0*rho.oldTime().internalField()
                       *vf.oldTime().internalField()*mesh().V0()
                      - coefft00*rho.oldTime().oldTime().internalField()
                       *vf.oldTime().oldTime().internalField()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*rho.boundaryField()*vf.boundaryField() -
                    (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
backwardDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT)*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0*vf.oldTime().internalField()*mesh().V0()
          - coefft00*vf.oldTime().oldTime().internalField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0*vf.oldTime().internalField()
          - coefft00*vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
backwardDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT*rho.value())*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*rho.value()*
        (
            coefft0*vf.oldTime().internalField()*mesh().V0()
          - coefft00*vf.oldTime().oldTime().internalField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*rho.value()*
        (
            coefft0*vf.oldTime().internalField()
          - coefft00*vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
backwardDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT)*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0*rho.oldTime().internalField()
           *vf.oldTime().internalField()*mesh().V0()
          - coefft00*rho.oldTime().oldTime().internalField()
           *vf.oldTime().oldTime().internalField()*mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0*rho.oldTime().internalField()
           *vf.oldTime().internalField()
          - coefft00*rho.oldTime().oldTime().internalField()
           *vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            rDeltaT*this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *(
                fvc::interpolate(rA)
               *(
                   coefft0*phi.oldTime()
                 - coefft00*phi.oldTime().oldTime()
                )
              - (
                    fvc::interpolate
                    (
                        rA*
                        (
                            coefft0*U.oldTime()
                          - coefft00*U.oldTime().oldTime()
                        )
                    ) & mesh().Sf()
                )
            )
        )
    );
}


template<class Type>
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ','
      + rho.name() + ','
      + U.name() + ','
      + phiAbs.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT*this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime())
               *(
                    coefft0*fvc::interpolate(rA*rho.oldTime())
                   *phiAbs.oldTime()
                  - coefft00*fvc::interpolate(rA*rho.oldTime().oldTime())
                   *phiAbs.oldTime().oldTime()
                  - (
                        fvc::interpolate
                        (
                            rA*
                            (
                                coefft0*rho.oldTime()*U.oldTime()
                              - coefft00*rho.oldTime().oldTime()
                               *U.oldTime().oldTime()
                            )
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                )
               *(
                    fvc::interpolate(rA*rho.oldTime())
                   *(
                       coefft0*phiAbs.oldTime()
                      /fvc::interpolate(rho.oldTime())
                     - coefft00*phiAbs.oldTime().oldTime()
                      /fvc::interpolate(rho.oldTime().oldTime())
                    )
                  - (
                        fvc::interpolate
                        (
                            rA*rho.oldTime()*
                            (
                                coefft0*U.oldTime()
                              - coefft00*U.oldTime().oldTime()
                            )
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phiAbs.oldTime())
               *(
                    fvc::interpolate(rA)
                   *(
                       coefft0*phiAbs.oldTime()
                     - coefft00*phiAbs.oldTime().oldTime()
                    )
                  - (
                        fvc::interpolate
                        (
                            rA*
                            (
                                coefft0*U.oldTime()
                              - coefft00*U.oldTime().oldTime()
                            )
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "backwardDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of phiAbs are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtConsistentPhiCorr
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(U);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    volScalarField V0oV
    (
        IOobject
        (
            "V0oV",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField V00oV
    (
        IOobject
        (
            "V00oV",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    );

    if (mesh().moving())
    {
        V0oV.internalField() = mesh().V0()/mesh().V();
        V0oV.correctBoundaryConditions();
    
        V00oV.internalField() = mesh().V00()/mesh().V();
        V00oV.correctBoundaryConditions();
    }

    // surfaceVectorField Uof = fvc::interpolate(U.oldTime()*V0oV);
    // surfaceVectorField Uoof = fvc::interpolate(U.oldTime().oldTime()*V00oV);

    fluxFieldType phiCorr
    (
        mesh().Sf()
      & (
            (
                coefft0*Uf.oldTime()*fvc::interpolate(V0oV)
              - coefft00*Uf.oldTime().oldTime()*fvc::interpolate(V00oV)
            )
          - (
                coefft0*fvc::interpolate(U.oldTime()*V0oV)
              - coefft00*fvc::interpolate(U.oldTime().oldTime()*V00oV)
            )
        )
    );

    forAll(phiCorr.boundaryField(), patchI)
    {
        if
        (
            U.boundaryField()[patchI].fixesValue()
         // || isA<slipFvPatchVectorField>(U.boundaryField()[patchI])
        )
        {
            phiCorr.boundaryField()[patchI] *= 0.0;
        }
    }

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            // this->fvcDdtPhiCoeff(U.oldTime(), (mesh().Sf() & Uf.oldTime()))
            rDeltaT*phiCorr*rAUf
        )
    );


    // const scalar deltaT = deltaT_();
    // const scalar deltaT0 = deltaT0_(U);

    // const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    // const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    // const scalar coefft0  = coefft + coefft00;

    // const scalar rDeltaT = 1.0/deltaT;

    // // Note: minus sign in gamma coefficient so we can simply add the fluxes
    // // together at the end
    // const dimensionedScalar beta("beta", dimless/dimTime, coefft0*rDeltaT);
    // const dimensionedScalar gamma("gamma", dimless/dimTime, -coefft00*rDeltaT);

    // // Calculate old and old-old flux contributions
    // fluxFieldType oldTimeFlux =
    //     beta*rAUf*(mesh().Sf() & faceU.oldTime());
    // fluxFieldType oldOldTimeFlux =
    //     gamma*rAUf*(mesh().Sf() & faceU.oldTime().oldTime());

    // if (mesh().moving())
    // {
    //     // Mesh is moving, need to take into account the ratio between old and
    //     // current cell volumes for old flux contribution
    //     volScalarField V0ByV
    //     (
    //         IOobject
    //         (
    //             "V0ByV",
    //             mesh().time().timeName(),
    //             mesh(),
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh(),
    //         dimensionedScalar("one", dimless, 1.0),
    //         zeroGradientFvPatchScalarField::typeName
    //     );
    //     V0ByV.internalField() = mesh().V0()/mesh().V();
    //     V0ByV.correctBoundaryConditions();

    //     // Correct old time flux contribution
    //     oldTimeFlux *= fvc::interpolate(V0ByV);


    //     // Also need to take into account the ratio between old-old and current
    //     // cell volumes for old-old time flux contribution
    //     volScalarField V00ByV
    //     (
    //         IOobject
    //         (
    //             "V00ByV",
    //             mesh().time().timeName(),
    //             mesh(),
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh(),
    //         dimensionedScalar("one", dimless, 1.0),
    //         zeroGradientFvPatchScalarField::typeName
    //     );
    //     V00ByV.internalField() = mesh().V00()/mesh().V();
    //     V00ByV.correctBoundaryConditions();

    //     // Correct old-old time flux contribution
    //     oldOldTimeFlux *= fvc::interpolate(V00ByV);
    // }

    // return oldTimeFlux + oldOldTimeFlux;
}


template<class Type>
tmp<surfaceScalarField> backwardDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    scalar coefftn_0 = 1 + coefft0_00;

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                mesh().phi().name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            coefftn_0*mesh().phi() - coefft0_00*mesh().phi().oldTime()
        )
    );

    // ZT: This bugfix is not correct. ESI version is correct!
    // Bugfix: missing possibility of having the variable time step
    // Reported by Sopheak Seng, Bureau Veritas, 6/Sep/2018.
    // const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    // const scalar coefft  = 1 + deltaT/(deltaT + deltaT0);

    // return coefft*mesh().phi() - coefft00*mesh().phi().oldTime();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
