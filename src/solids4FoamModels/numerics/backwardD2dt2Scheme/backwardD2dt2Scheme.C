/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "backwardD2dt2Scheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "backwardDdtScheme.H"
#include "EulerD2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT_() const
{
    return mesh().time().deltaT().value();
}


template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0().value();
}


template<class Type>
scalar backwardD2dt2Scheme<Type>::deltaT0_
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if
    (
        vf.oldTime().oldTime().timeIndex()
     == vf.oldTime().oldTime().oldTime().timeIndex()
    )
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
backwardD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first timne step
    if (mesh().time().timeIndex() == 1)
    {
        return EulerD2dt2Scheme<Type>(mesh()).fvcD2dt2(vf);
    }

    IOobject d2dt2IOobject
    (
        "d2dt2(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar coefft = 1.5;
    const scalar coefft0 = 2.0;
    const scalar coefft00 = 0.5;

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            d2dt2IOobject,
            rDeltaT*
            (
                coefft*backwardDdtScheme<Type>
                (
                    mesh()
                ).fvcDdt(vf)
              - coefft0*backwardDdtScheme<Type>
                (
                    mesh()
                ).fvcDdt(vf.oldTime())
              + coefft00*backwardDdtScheme<Type>
                (
                    mesh()
                ).fvcDdt(vf.oldTime().oldTime())
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
backwardD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first timne step
    if (mesh().time().timeIndex() == 1)
    {
        return EulerD2dt2Scheme<Type>(mesh()).fvcD2dt2(rho, vf);
    }

    IOobject d2dt2IOobject
    (
        "d2dt2(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar coefft = 1.5;
    const scalar coefft0 = 2.0;
    const scalar coefft00 = 0.5;

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            d2dt2IOobject,
            rDeltaT*
            (
                coefft*rho*backwardDdtScheme<Type>
                (
                    mesh()
                ).fvcDdt(vf)
                - coefft0*rho.oldTime()*backwardDdtScheme<Type>
                (
                    mesh()
                ).fvcDdt(vf.oldTime())
                + coefft00*rho.oldTime().oldTime()*backwardDdtScheme<Type>
                (
                    mesh()
                ).fvcDdt(vf.oldTime().oldTime())
            )
        )
    );
}


template<class Type>
tmp<fvMatrix<Type> >
backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first timne step
    if (mesh().time().timeIndex() == 1)
    {
        return EulerD2dt2Scheme<Type>(mesh()).fvmD2dt2(vf);
    }

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    const scalar rDeltaT = 1.0/deltaT_();
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    fvm = coefft*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
       *backwardDdtScheme<Type>(mesh()).fvmDdt(vf);

#ifdef FOAMEXTEND
    fvm.source() += rDeltaT*mesh().V()*
    (
        coefft0*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().internalField()
      - coefft00*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().internalField()
    );
#else
    fvm.source() += rDeltaT*mesh().V()*
    (
        coefft0*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().primitiveField()
      - coefft00*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().primitiveField()
    );
#endif

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first timne step
    if (mesh().time().timeIndex() == 1)
    {
        return EulerD2dt2Scheme<Type>(mesh()).fvmD2dt2(rho, vf);
    }

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    const scalar rDeltaT = 1.0/deltaT_();
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    fvm = coefft*rho*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
       *backwardDdtScheme<Type>(mesh()).fvmDdt(vf);

#ifdef FOAMEXTEND
    fvm.source() += rDeltaT*rho*mesh().V()*
    (
        coefft0*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().internalField()
      - coefft00*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().internalField()
    );
#else
    fvm.source() += rDeltaT*rho*mesh().V()*
    (
        coefft0*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().primitiveField()
      - coefft00*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().primitiveField()
    );
#endif

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
backwardD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (mag(mesh().time().deltaT() - mesh().time().deltaT0()).value() > SMALL)
    {
        notImplemented
        (
            "backwardD2dt2Scheme not implemented for variable time steps"
        );
    }

    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    // Default to 1st order Euler on the first timne step
    if (mesh().time().timeIndex() == 1)
    {
        return EulerD2dt2Scheme<Type>(mesh()).fvmD2dt2(rho, vf);
    }

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    const scalar rDeltaT = 1.0/deltaT_();
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    fvm = coefft*rho*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
       *backwardDdtScheme<Type>(mesh()).fvmDdt(vf);

#ifdef FOAMEXTEND
    fvm.source() += rDeltaT*mesh().V()*
    (
        coefft0*rho.oldTime()*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().internalField()
      - coefft00*rho.oldTime().oldTime()*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().internalField()
    );
#else
    fvm.source() += rDeltaT*mesh().V()*
    (
        coefft0*rho.oldTime()*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime())().primitiveField()
      - coefft00*rho.oldTime().oldTime()*backwardDdtScheme<Type>
        (
            mesh()
        ).fvcDdt(vf.oldTime().oldTime())().primitiveField()
    );
#endif

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
