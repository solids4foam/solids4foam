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

Description
    Implementation of the BDF ddt scheme: the coefficients of each order,
    the start-up handling, the fvm::ddt and fvc::ddt members and the ddtCorr
    flux corrections.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "bdfDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
label bdfDdtScheme<Type>::readOrder(Istream& is)
{
    // Check for the end of the stream before reading: an ITstream reports
    // eof as soon as its last token has been read
    if (is.eof())
    {
        FatalIOErrorInFunction(is)
            << "The BDF ddt scheme requires an order (1 to 6): "
            << "e.g. \"BDF 2;\" in ddtSchemes"
            << exit(FatalIOError);
    }

    token t(is);

    if (!t.isLabel())
    {
        FatalIOErrorInFunction(is)
            << "The BDF ddt scheme requires an integer order (1 to 6), "
            << "but found: " << t.info()
            << exit(FatalIOError);
    }

    const label order = t.labelToken();

    if (order < 1 || order > 6)
    {
        FatalIOErrorInFunction(is)
            << "The order of the BDF ddt scheme must be 1 to 6, "
            << "but it is " << order << exit(FatalIOError);
    }

    return order;
}


template<class Type>
bdfDdtScheme<Type>::bdfDdtScheme(const fvMesh& mesh, Istream& is)
:
    ddtScheme<Type>(mesh, is),
    order_(readOrder(is))
{
    if (is.good() && !is.eof())
    {
        is >> this->ddtPhiCoeff_;
    }

    if (order_ >= 3 && !bdfDdtSchemeWarned)
    {
        bdfDdtSchemeWarned = true;

        WarningInFunction
            << "The " << type() << " " << order_ << " ddt scheme is not "
            << "A-stable, only A(alpha)-stable, and its stability region is "
            << "smaller for higher orders. It suits diffusion-dominated, "
            << "stiff problems; in strongly convection-dominated flow at "
            << "large Courant numbers it can be less robust than backward, "
            << "and BDF 5 and 6 are not recommended for it"
            << endl;
    }

    checkStaticMesh();
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

template<class Type>
void bdfDdtScheme<Type>::checkStaticMesh() const
{
    if (mesh().moving())
    {
        FatalErrorInFunction
            << "BDF " << order_ << " ddt is implemented for static meshes"
            << " only: fvMesh stores V0 and V00, fewer old cell volumes than"
            << " an order-" << order_ << " formula needs."
            << exit(FatalError);
    }
}


template<class Type>
template<class GeoField>
const GeoField& bdfDdtScheme<Type>::level(const GeoField& vf, const label k)
{
    const GeoField* f = &vf;

    for (label i = 0; i < k; ++i)
    {
        f = &f->oldTime();
    }

    return *f;
}


template<class Type>
template<class GeoField>
void bdfDdtScheme<Type>::touch(const GeoField& vf) const
{
    level(vf, order_);
}


template<class Type>
template<class GeoField>
label bdfDdtScheme<Type>::order(const GeoField& vf) const
{
    // Create the whole chain first. Accessing level 1 also shifts the old
    // times when this is a new time step, so the indices compared below are
    // those of the step being solved.
    touch(vf);

    label q = 1;

    while
    (
        q < order_
     && level(vf, q + 1).timeIndex() != level(vf, q).timeIndex()
    )
    {
        ++q;
    }

    return q;
}


template<class Type>
template<class GeoField1, class GeoField2>
label bdfDdtScheme<Type>::order(const GeoField1& f1, const GeoField2& f2) const
{
    return min(order(f1), order(f2));
}


template<class Type>
template<class GeoField1, class GeoField2, class GeoField3>
label bdfDdtScheme<Type>::order
(
    const GeoField1& f1,
    const GeoField2& f2,
    const GeoField3& f3
) const
{
    return min(order(f1, f2), order(f3));
}


template<class Type>
scalar bdfDdtScheme<Type>::deltaT0() const
{
    const Time& runTime = mesh().time();
    const scalar deltaT0 = runTime.deltaT0Value();

    // On the first step after a restart with adjustTimeStep no, Time has set
    // deltaT0 to the deltaT of controlDict, not the step the old levels were
    // stored with. That step is in <startTime>/uniform/time, which the
    // previous run wrote.
    if (runTime.timeIndex() == runTime.startTimeIndex() + 1)
    {
        const IOdictionary timeDict
        (
            IOobject
            (
                "time",
                runTime.timeName(runTime.startTime().value()),
                "uniform",
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        scalar storedDeltaT = deltaT0;

        // The stored value has only writePrecision digits: use it only if
        // it differs by more than that, so that an unchanged step keeps the
        // exact deltaT0 and a restart reproduces a continuous run
        if
        (
            timeDict.readIfPresent("deltaT", storedDeltaT)
         && mag(storedDeltaT - deltaT0) > deltaTTol_*deltaT0
        )
        {
            return storedDeltaT;
        }
    }

    return deltaT0;
}


template<class Type>
scalarList bdfDdtScheme<Type>::weights(const label q) const
{
    if (q > 1)
    {
        const scalar deltaT = mesh().time().deltaTValue();
        const scalar deltaT0 = this->deltaT0();

        if (order_ == 2)
        {
            // Variable-step BDF2, which needs only deltaT and deltaT0. These
            // are backward's coefficients, and reduce to 3/2, -2, 1/2 for a
            // constant step.
            const scalar omega = deltaT/deltaT0;

            return scalarList
            (
                {
                    (1 + 2*omega)/(1 + omega),
                    -(1 + omega),
                    sqr(omega)/(1 + omega)
                }
            );
        }

        // From order 3 the step ratios of older steps are needed as well,
        // which Time does not record, so the step must be constant. This
        // also applies to the BDF2 start-up steps of a higher order, whose
        // levels are reused by the full-order formula later.
        if (mag(deltaT - deltaT0) > deltaTTol_*deltaT)
        {
            FatalErrorInFunction
                << "BDF " << order_ << " ddt requires a constant time step,"
                << " but deltaT changed from " << deltaT0 << " to " << deltaT
                << " at time " << mesh().time().timeName() << "." << nl
                << "    Its coefficients are for a uniform step; applying"
                << " them to a non-uniform stencil would silently lose the"
                << " order of accuracy. BDF 1 and BDF 2 support a variable"
                << " time step, including a change of deltaT on restart."
                << exit(FatalError);
        }
    }

    switch (q)
    {
        case 1:
            return scalarList({1, -1});

        case 2:
            return scalarList({3.0/2.0, -2, 1.0/2.0});

        case 3:
            return scalarList({11.0/6.0, -3, 3.0/2.0, -1.0/3.0});

        case 4:
            return scalarList({25.0/12.0, -4, 3, -4.0/3.0, 1.0/4.0});

        case 5:
            return scalarList
            (
                {137.0/60.0, -5, 5, -10.0/3.0, 5.0/4.0, -1.0/5.0}
            );

        case 6:
            return scalarList
            (
                {
                    49.0/20.0, -6, 15.0/2.0, -20.0/3.0, 15.0/4.0, -6.0/5.0,
                    1.0/6.0
                }
            );

        default:
            FatalErrorInFunction
                << "No BDF coefficients for order " << q
                << exit(FatalError);
    }

    return scalarList();
}


template<class Type>
template<class GeoField>
tmp<GeoField> bdfDdtScheme<Type>::history
(
    const GeoField& vf,
    const scalarList& a
) const
{
    touch(vf);

    tmp<GeoField> tsum(-a[1]*level(vf, 1));

    for (label k = 2; k < a.size(); ++k)
    {
        tsum.ref() -= a[k]*level(vf, k);
    }

    return tsum;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> bdfDdtScheme<Type>::history
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const scalarList& a
) const
{
    touch(rho);
    touch(vf);

    tmp<GeometricField<Type, fvPatchField, volMesh>> tsum
    (
        -a[1]*level(rho, 1)*level(vf, 1)
    );

    for (label k = 2; k < a.size(); ++k)
    {
        tsum.ref() -= a[k]*level(rho, k)*level(vf, k);
    }

    return tsum;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> bdfDdtScheme<Type>::history
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const scalarList& a
) const
{
    touch(alpha);
    touch(rho);
    touch(vf);

    tmp<GeometricField<Type, fvPatchField, volMesh>> tsum
    (
        -a[1]*level(alpha, 1)*level(rho, 1)*level(vf, 1)
    );

    for (label k = 2; k < a.size(); ++k)
    {
        tsum.ref() -= a[k]*level(alpha, k)*level(rho, k)*level(vf, k);
    }

    return tsum;
}


template<class Type>
template<class GeoField, class Reconstruct>
tmp<GeoField> bdfDdtScheme<Type>::fluxHistory
(
    const GeoField& flux,
    const scalarList& a,
    const Reconstruct& reconstruct
) const
{
    // Levels 1 to stored are genuine; older ones are copies of the level
    // above, made because the flux was not stored that far back
    const label stored = order(flux);

    tmp<GeoField> tsum(-a[1]*level(flux, 1));

    for (label k = 2; k < a.size(); ++k)
    {
        if (k <= stored)
        {
            tsum.ref() -= a[k]*level(flux, k);
        }
        else
        {
            tsum.ref() -= a[k]*reconstruct(k);
        }
    }

    return tsum;
}


template<class Type>
IOobject bdfDdtScheme<Type>::ddtIOobject(const word& name) const
{
    return IOobject
    (
        "ddt(" + name + ')',
        mesh().time().timeName(),
        mesh()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdfDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    // The mesh is static, so the time derivative of a constant is zero.
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        ddtIOobject(dt.name()),
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero),
        calculatedFvPatchField<Type>::typeName
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdfDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(vf)));

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject(vf.name()),
            rDeltaT*(a[0]*vf - history(vf, a))
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdfDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(vf)));

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject(rho.name() + ',' + vf.name()),
            rDeltaT*rho*(a[0]*vf - history(vf, a))
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdfDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(rho, vf)));

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject(rho.name() + ',' + vf.name()),
            rDeltaT*(a[0]*rho*vf - history(rho, vf, a))
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdfDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(alpha, rho, vf)));

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject(alpha.name() + ',' + rho.name() + ',' + vf.name()),
            rDeltaT*(a[0]*alpha*rho*vf - history(alpha, rho, vf, a))
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
bdfDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();
    const scalarList a(weights(order(vf)));

    fvm.diag() = (a[0]*rDeltaT)*mesh().V();

    fvm.source() = rDeltaT*mesh().V()*history(vf, a)().primitiveField();

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
bdfDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();
    const scalarList a(weights(order(vf)));

    fvm.diag() = (a[0]*rDeltaT*rho.value())*mesh().V();

    fvm.source() =
        rDeltaT*rho.value()*mesh().V()*history(vf, a)().primitiveField();

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
bdfDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();
    const scalarList a(weights(order(rho, vf)));

    fvm.diag() = (a[0]*rDeltaT)*rho.primitiveField()*mesh().V();

    fvm.source() =
        rDeltaT*mesh().V()*history(rho, vf, a)().primitiveField();

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
bdfDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();
    const scalarList a(weights(order(alpha, rho, vf)));

    fvm.diag() =
        (a[0]*rDeltaT)*alpha.primitiveField()*rho.primitiveField()*mesh().V();

    fvm.source() =
        rDeltaT*mesh().V()*history(alpha, rho, vf, a)().primitiveField();

    return tfvm;
}


template<class Type>
tmp<typename bdfDdtScheme<Type>::fluxFieldType>
bdfDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(U)));

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
            this->fvcDdtPhiCoeff(U.oldTime(), (mesh().Sf() & Uf.oldTime()))
           *rDeltaT
           *(
                mesh().Sf()
              & (
                    fluxHistory
                    (
                        Uf,
                        a,
                        [&](const label k)
                        {
                            return fvc::interpolate(level(U, k));
                        }
                    )
                  - fvc::interpolate(history(U, a))
                )
            )
        )
    );
}


template<class Type>
tmp<typename bdfDdtScheme<Type>::fluxFieldType>
bdfDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(U)));

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *rDeltaT
           *(
                fluxHistory
                (
                    phi,
                    a,
                    [&](const label k)
                    {
                        return fvc::dotInterpolate(mesh().Sf(), level(U, k));
                    }
                )
              - fvc::dotInterpolate(mesh().Sf(), history(U, a))
            )
        )
    );
}


template<class Type>
tmp<typename bdfDdtScheme<Type>::fluxFieldType>
bdfDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(rho, U)));

    const word name
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')'
    );

    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        // U is a velocity and Uf a momentum: form the old momenta rho*U
        const GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject(name, mesh().time().timeName(), mesh()),
                this->fvcDdtPhiCoeff
                (
                    rhoU0,
                    mesh().Sf() & Uf.oldTime(),
                    rho.oldTime()
                )
               *rDeltaT
               *(
                    mesh().Sf()
                  & (
                        fluxHistory
                        (
                            Uf,
                            a,
                            [&](const label k)
                            {
                                return fvc::interpolate
                                (
                                    level(rho, k)*level(U, k)
                                );
                            }
                        )
                      - fvc::interpolate(history(rho, U, a))
                    )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject(name, mesh().time().timeName(), mesh()),
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    mesh().Sf() & Uf.oldTime(),
                    rho.oldTime()
                )
               *rDeltaT
               *(
                    mesh().Sf()
                  & (
                        fluxHistory
                        (
                            Uf,
                            a,
                            [&](const label k)
                            {
                                return fvc::interpolate(level(U, k));
                            }
                        )
                      - fvc::interpolate(history(U, a))
                    )
                )
            )
        );
    }

    FatalErrorInFunction
        << "dimensions of Uf are not correct"
        << abort(FatalError);

    return fluxFieldType::null();
}


template<class Type>
tmp<typename bdfDdtScheme<Type>::fluxFieldType>
bdfDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const scalarList a(weights(order(rho, U)));

    const word name
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')'
    );

    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        // U is a velocity and phi a mass flux: form the old momenta rho*U
        const GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject(name, mesh().time().timeName(), mesh()),
                this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), rho.oldTime())
               *rDeltaT
               *(
                    fluxHistory
                    (
                        phi,
                        a,
                        [&](const label k)
                        {
                            return fvc::dotInterpolate
                            (
                                mesh().Sf(),
                                level(rho, k)*level(U, k)
                            );
                        }
                    )
                  - fvc::dotInterpolate(mesh().Sf(), history(rho, U, a))
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject(name, mesh().time().timeName(), mesh()),
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), rho.oldTime())
               *rDeltaT
               *(
                    fluxHistory
                    (
                        phi,
                        a,
                        [&](const label k)
                        {
                            return fvc::dotInterpolate
                            (
                                mesh().Sf(),
                                level(U, k)
                            );
                        }
                    )
                  - fvc::dotInterpolate(mesh().Sf(), history(U, a))
                )
            )
        );
    }

    FatalErrorInFunction
        << "dimensions of phi are not correct"
        << abort(FatalError);

    return fluxFieldType::null();
}


template<class Type>
tmp<surfaceScalarField> bdfDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    // Unreachable in practice: construction already fails on a moving mesh,
    // the only case in which the mesh flux is asked for. The name is given
    // literally, as mesh().phi() is fatal on a static mesh.
    checkStaticMesh();

    return surfaceScalarField::New
    (
        "meshPhi",
        IOobject::NO_REGISTER,
        mesh(),
        dimensionedScalar(dimVolume/dimTime, Zero)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
