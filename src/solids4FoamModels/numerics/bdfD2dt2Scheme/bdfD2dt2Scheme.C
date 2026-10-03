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
tmp<Field<Type> > bdfD2dt2Scheme<Type>::levelSum
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const label j,
    const label i0,
    const label p
) const
{
    const scalarList a = coeffs(p);

    tmp<Field<Type> > tsum(new Field<Type>(vf.size(), pTraits<Type>::zero));
#ifdef FOAMEXTEND
    Field<Type>& sum = tsum();
#else
    Field<Type>& sum = tsum.ref();
#endif

    for (label i = i0; i <= p; i++)
    {
#ifdef FOAMEXTEND
        sum += a[i]*oldTimeLevel(vf, j + i).internalField();
#else
        sum += a[i]*oldTimeLevel(vf, j + i).primitiveField();
#endif
    }

    return tsum;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
bdfD2dt2Scheme<Type>::levelDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const label j,
    const label p
) const
{
    const scalarList a = coeffs(p);
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
    const GeometricField<Type, fvPatchField, volMesh>& vfj =
        oldTimeLevel(vf, j);

    tmp<GeometricField<Type, fvPatchField, volMesh> > tddt
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "ddt(" + vfj.name() + ')',
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rDeltaT*a[0]*vfj
        )
    );
#ifdef FOAMEXTEND
    GeometricField<Type, fvPatchField, volMesh>& ddt = tddt();
#else
    GeometricField<Type, fvPatchField, volMesh>& ddt = tddt.ref();
#endif

    for (label i = 1; i <= p; i++)
    {
        ddt += rDeltaT*a[i]*oldTimeLevel(vf, j + i);
    }

    return tddt;
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
    order_(readLabel(is))
{
    if (order_ < 1 || order_ > 6)
    {
        FatalIOErrorInFunction(is)
            << "The order of the BDF d2dt2 scheme must be 1 to 6, "
            << "but it is " << order_ << exit(FatalIOError);
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
    const label p = levelOrder(n, 0);
    const scalarList a = coeffs(p);
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

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
            rDeltaT*a[0]*levelDdt(vf, 0, levelOrder(n, 0))
        )
    );
#ifdef FOAMEXTEND
    GeometricField<Type, fvPatchField, volMesh>& d2dt2 = td2dt2();
#else
    GeometricField<Type, fvPatchField, volMesh>& d2dt2 = td2dt2.ref();
#endif

    for (label j = 1; j <= p; j++)
    {
        d2dt2 += rDeltaT*a[j]*levelDdt(vf, j, levelOrder(n, j));
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
    const label p = levelOrder(n, 0);
    const scalarList a = coeffs(p);
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

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
            rDeltaT*a[0]*rho*levelDdt(vf, 0, levelOrder(n, 0))
        )
    );
#ifdef FOAMEXTEND
    GeometricField<Type, fvPatchField, volMesh>& d2dt2 = td2dt2();
#else
    GeometricField<Type, fvPatchField, volMesh>& d2dt2 = td2dt2.ref();
#endif

    for (label j = 1; j <= p; j++)
    {
        d2dt2 += rDeltaT*a[j]*oldTimeLevel(rho, j)*levelDdt(vf, j, levelOrder(n, j));
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

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    const label n = nOldTimes(vf);
    const label p = levelOrder(n, 0);
    const scalarList a = coeffs(p);
    const scalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT().value());
    const scalarField& V = mesh().V();

    // Old-time terms: the current velocity uses levels 1 to p and the
    // velocity at old-time level j uses levels j to j + levelOrder(n, j)
    Field<Type> oldTerms(a[0]*levelSum(vf, 0, 1, levelOrder(n, 0)));

    for (label j = 1; j <= p; j++)
    {
        oldTerms += a[j]*levelSum(vf, j, 0, levelOrder(n, j));
    }

    fvm.diag() = sqr(a[0])*rDeltaT2*V;
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
    checkMesh();

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

    const label n = nOldTimes(vf);
    const label p = levelOrder(n, 0);
    const scalarList a = coeffs(p);
    const scalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT().value());
    const scalarField& V = mesh().V();

    Field<Type> oldTerms(a[0]*levelSum(vf, 0, 1, levelOrder(n, 0)));

    for (label j = 1; j <= p; j++)
    {
        oldTerms += a[j]*levelSum(vf, j, 0, levelOrder(n, j));
    }

    fvm.diag() = sqr(a[0])*rDeltaT2*rho.value()*V;
    fvm.source() = -rDeltaT2*rho.value()*V*oldTerms;

    return tfvm;
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

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    const label n = nOldTimes(vf);
    const label p = levelOrder(n, 0);
    const scalarList a = coeffs(p);
    const scalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT().value());
    const scalarField& V = mesh().V();

    // Create all density old-time levels, so that they are stored from the
    // first time-step
    oldTimeLevel(rho, order_);

    // The density at old-time level j multiplies the velocity at level j
#ifdef FOAMEXTEND
    const scalarField& rhoI = rho.internalField();
#else
    const scalarField& rhoI = rho.primitiveField();
#endif

    Field<Type> oldTerms(a[0]*rhoI*levelSum(vf, 0, 1, levelOrder(n, 0)));

    for (label j = 1; j <= p; j++)
    {
#ifdef FOAMEXTEND
        const scalarField& rhoJ = oldTimeLevel(rho, j).internalField();
#else
        const scalarField& rhoJ = oldTimeLevel(rho, j).primitiveField();
#endif

        oldTerms += a[j]*rhoJ*levelSum(vf, j, 0, levelOrder(n, j));
    }

    fvm.diag() = sqr(a[0])*rDeltaT2*rhoI*V;
    fvm.source() = -rDeltaT2*V*oldTerms;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
