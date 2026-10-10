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

\*---------------------------------------------------------------------------*/

#include "bdfD2dt2Scheme.H"
#include "fvMatrices.H"
#include "token.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalarList bdfD2dt2Scheme<Type>::coeffs(const label order)
{
    scalarList a(order + 1);

    if (order == 1)
    {
        a[0] = 1.0;
        a[1] = -1.0;
    }
    else if (order == 2)
    {
        a[0] = 3.0/2.0;
        a[1] = -2.0;
        a[2] = 1.0/2.0;
    }
    else if (order == 3)
    {
        a[0] = 11.0/6.0;
        a[1] = -3.0;
        a[2] = 3.0/2.0;
        a[3] = -1.0/3.0;
    }
    else if (order == 4)
    {
        a[0] = 25.0/12.0;
        a[1] = -4.0;
        a[2] = 3.0;
        a[3] = -4.0/3.0;
        a[4] = 1.0/4.0;
    }
    else if (order == 5)
    {
        a[0] = 137.0/60.0;
        a[1] = -5.0;
        a[2] = 5.0;
        a[3] = -10.0/3.0;
        a[4] = 5.0/4.0;
        a[5] = -1.0/5.0;
    }
    else if (order == 6)
    {
        a[0] = 147.0/60.0;
        a[1] = -6.0;
        a[2] = 15.0/2.0;
        a[3] = -20.0/3.0;
        a[4] = 15.0/4.0;
        a[5] = -6.0/5.0;
        a[6] = 1.0/6.0;
    }
    else
    {
        FatalErrorInFunction
            << "Order " << order << " is not implemented"
            << abort(FatalError);
    }

    return a;
}


template<class Type>
template<class GeoField>
const GeoField& bdfD2dt2Scheme<Type>::oldTimeLevel
(
    const GeoField& vf,
    const label n
)
{
    const GeoField* fieldPtr = &vf;

    for (label i = 0; i < n; i++)
    {
        fieldPtr = &fieldPtr->oldTime();
    }

    return *fieldPtr;
}


template<class Type>
label bdfD2dt2Scheme<Type>::nOldTimes
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Count the genuine old-time levels. A newly created old-time field is a
    // copy of its parent with the same time index. All 2*order_ levels are
    // accessed so that they are created, and stored, from the first
    // time-step
    label nOldTimes = 0;
    const GeometricField<Type, fvPatchField, volMesh>* fieldPtr = &vf;

    for (label i = 1; i <= 2*order_; i++)
    {
        const GeometricField<Type, fvPatchField, volMesh>& field0 =
            fieldPtr->oldTime();

        if (nOldTimes == i - 1 && field0.timeIndex() != fieldPtr->timeIndex())
        {
            nOldTimes = i;
        }

        fieldPtr = &field0;
    }

    // The coefficients assume a constant time-step
    if (nOldTimes > 1)
    {
        const scalar deltaT = mesh().time().deltaT().value();
        const scalar deltaT0 = mesh().time().deltaT0().value();

        if (mag(deltaT - deltaT0) > 1e-5*deltaT)
        {
            FatalErrorInFunction
                << "The " << type() << " d2dt2 scheme requires a constant "
                << "time-step, but the time-step changed from " << deltaT0
                << " to " << deltaT << nl
                << "Use the Euler or backward d2dt2 scheme for a variable "
                << "time-step" << exit(FatalError);
        }
    }

    return nOldTimes;
}


template<class Type>
label bdfD2dt2Scheme<Type>::levelOrder
(
    const label nOldTimes,
    const label j
) const
{
    // The velocity at old-time level j was the current velocity when
    // nOldTimes - j genuine old-time levels were available, so it uses the
    // order it had then. A velocity with no genuine old-time levels, at the
    // initial time, is zero, as in the backward scheme
    return max(1, min(order_, nOldTimes - j));
}


template<class Type>
scalarList bdfD2dt2Scheme<Type>::composedCoeffs
(
    const label nOldTimes
) const
{
    const label p0 = levelOrder(nOldTimes, 0);
    const scalarList a0 = coeffs(p0);

    scalarList C(2*order_ + 1, 0.0);

    for (label j = 0; j <= p0; j++)
    {
        const label pj = levelOrder(nOldTimes, j);
        const scalarList aj = coeffs(pj);

        for (label i = 0; i <= pj; i++)
        {
            C[j + i] += a0[j]*aj[i];
        }
    }

    return C;
}


template<class Type>
label bdfD2dt2Scheme<Type>::readOrder(Istream& is)
{
    // Check for the end of the stream before reading: an ITstream reports
    // eof as soon as its last token has been read
    if (is.eof())
    {
        FatalIOErrorInFunction(is)
            << "The BDF d2dt2 scheme requires an order (1 to 6): "
            << "e.g. \"BDF 2;\" in d2dt2Schemes"
            << exit(FatalIOError);
    }

    token t(is);

    if (!t.isLabel())
    {
        FatalIOErrorInFunction(is)
            << "The BDF d2dt2 scheme requires an integer order (1 to 6), "
            << "but found: " << t.info()
            << exit(FatalIOError);
    }

    const label order = t.labelToken();

    if (order < 1 || order > 6)
    {
        FatalIOErrorInFunction(is)
            << "The order of the BDF d2dt2 scheme must be 1 to 6, "
            << "but it is " << order << exit(FatalIOError);
    }

    return order;
}


template<class Type>
void bdfD2dt2Scheme<Type>::checkMesh() const
{
    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
bdfD2dt2Scheme<Type>::bdfD2dt2Scheme(const fvMesh& mesh, Istream& is)
:
    d2dt2Scheme<Type>(mesh, is),
    order_(readOrder(is))
{
    static bool warned = false;

    if (order_ >= 3 && !warned)
    {
        warned = true;

        WarningInFunction
            << "The " << type() << " " << order_ << " d2dt2 scheme "
            << "adds energy to undamped modes for finite time steps"
            << (order_ >= 5 ? " and can diverge on under-resolved modes" : "")
            << "; it should be used only where the problem provides "
            << "damping. For undamped problems, use backward or NewmarkBeta"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
bdfD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    checkMesh();

    const label n = nOldTimes(vf);
    const scalarList C = composedCoeffs(n);
    const dimensionedScalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT());

    tmp<GeometricField<Type, fvPatchField, volMesh> > td2dt2
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "d2dt2(" + vf.name() + ')',
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rDeltaT2*C[0]*vf
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& d2dt2 = tmpRef(td2dt2);

    const GeometricField<Type, fvPatchField, volMesh>* fieldPtr = &vf;

    for (label m = 1; m <= 2*order_; m++)
    {
        fieldPtr = &fieldPtr->oldTime();

        if (mag(C[m]) > 1e-12)
        {
            d2dt2 += rDeltaT2*C[m]*(*fieldPtr);
        }
    }

    return td2dt2;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
bdfD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    checkMesh();

    const label n = nOldTimes(vf);
    const label p0 = levelOrder(n, 0);
    const scalarList a0 = coeffs(p0);
    const dimensionedScalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT());

    // Create all density old-time levels, so that they are stored from the
    // first time-step
    oldTimeLevel(rho, order_);

    tmp<GeometricField<Type, fvPatchField, volMesh> > td2dt2
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "d2dt2(" + rho.name() + ',' + vf.name() + ')',
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rDeltaT2*sqr(a0[0])*rho*vf
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& d2dt2 = tmpRef(td2dt2);

    for (label i = 1; i <= p0; i++)
    {
        d2dt2 += rDeltaT2*a0[0]*a0[i]*rho*oldTimeLevel(vf, i);
    }

    for (label j = 1; j <= p0; j++)
    {
        const label pj = levelOrder(n, j);
        const scalarList aj = coeffs(pj);
        const volScalarField& rhoJ = oldTimeLevel(rho, j);

        for (label i = 0; i <= pj; i++)
        {
            d2dt2 += rDeltaT2*a0[j]*aj[i]*rhoJ*oldTimeLevel(vf, j + i);
        }
    }

    return td2dt2;
}


template<class Type>
tmp<fvMatrix<Type> >
bdfD2dt2Scheme<Type>::fvmD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    checkMesh();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tmpRef(tfvm);

    const label n = nOldTimes(vf);
    const scalarList C = composedCoeffs(n);
    const scalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT().value());
    const scalarField& V = mesh().V();

    Field<Type> oldTerms(vf.size(), pTraits<Type>::zero);
    const GeometricField<Type, fvPatchField, volMesh>* fieldPtr = &vf;

    for (label m = 1; m <= 2*order_; m++)
    {
        fieldPtr = &fieldPtr->oldTime();

        if (mag(C[m]) > 1e-12)
        {
            oldTerms += C[m]*primitiveField(*fieldPtr);
        }
    }

    fvm.diag() = C[0]*rDeltaT2*V;
    fvm.source() = -rDeltaT2*V*oldTerms;

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
bdfD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return rho*fvmD2dt2(vf);
}


template<class Type>
tmp<fvMatrix<Type> >
bdfD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    checkMesh();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tmpRef(tfvm);

    const label n = nOldTimes(vf);
    const label p0 = levelOrder(n, 0);
    const scalarList a0 = coeffs(p0);
    const scalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT().value());
    const scalarField& V = mesh().V();

    // Create all density old-time levels, so that they are stored from the
    // first time-step
    oldTimeLevel(rho, order_);

    const scalarField& rhoI = primitiveField(rho);

    Field<Type> oldTerms(vf.size(), pTraits<Type>::zero);

    for (label i = 1; i <= p0; i++)
    {
        const Field<Type>& vfi = primitiveField(oldTimeLevel(vf, i));
        oldTerms += a0[0]*a0[i]*rhoI*vfi;
    }

    for (label j = 1; j <= p0; j++)
    {
        const label pj = levelOrder(n, j);
        const scalarList aj = coeffs(pj);
        const scalarField& rhoJ = primitiveField(oldTimeLevel(rho, j));

        for (label i = 0; i <= pj; i++)
        {
            const Field<Type>& vfji = primitiveField(oldTimeLevel(vf, j + i));
            oldTerms += a0[j]*aj[i]*rhoJ*vfji;
        }
    }

    fvm.diag() = sqr(a0[0])*rDeltaT2*rhoI*V;
    fvm.source() = -rDeltaT2*V*oldTerms;

    return tfvm;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
