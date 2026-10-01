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

#include "NewmarkBetaD2dt2Scheme.H"
#include "fvMatrices.H"
#include "calculatedFvPatchFields.H"
#include "compatibilityFunctions.H"
#include "fvcD2dt2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
word NewmarkBetaD2dt2Scheme<Type>::startTimeName() const
{
    // The name of the start time directory, rather than the start time
    // written with the current timePrecision, which differs from it if the
    // directory was written with more digits
    const Time& runTime = mesh().time();
    const scalar startTime = runTime.startTime().value();
    const instant closest = runTime.findClosestTime(startTime);

    if (mag(closest.value() - startTime) < 0.5*runTime.deltaTValue())
    {
        return closest.name();
    }

    return runTime.timeName(startTime);
}


template<class Type>
GeometricField<Type, fvPatchField, volMesh>&
NewmarkBetaD2dt2Scheme<Type>::stateField
(
    const word& name,
    const dimensioned<Type>& value
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& mesh = this->mesh();

    if (mesh.foundObject<fieldType>(name))
    {
        return const_cast<fieldType&>(mesh.lookupObject<fieldType>(name));
    }

    IOobject io
    (
        name,
        startTimeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

#ifdef OPENFOAM_NOT_EXTEND
    const bool present = io.typeHeaderOk<fieldType>(true);
#else
    const bool present = io.headerOk();
#endif

    if (present)
    {
        Info<< type() << ": reading " << name << endl;

        fieldType* fieldPtr = new fieldType(io, mesh);

        return regIOobject::store(fieldPtr);
    }

    io.readOpt() = IOobject::NO_READ;

    return regIOobject::store
    (
        new fieldType
        (
            io,
            mesh,
            value,
            calculatedFvPatchField<Type>::typeName
        )
    );
}


template<class Type>
void NewmarkBetaD2dt2Scheme<Type>::updateState
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (mesh().moving())
    {
        notImplemented(type() + ": not implemented for a moving mesh");
    }

    const word& vfName = vf.name();

    if (vfName.size() > 2 && vfName.substr(vfName.size() - 2) == "_0")
    {
        FatalErrorInFunction
            << "The " << type() << " d2dt2 scheme cannot be applied to the "
            << "old-time field " << vfName << ", because its velocity and "
            << "acceleration are stored for the current field only." << nl
            << "The updated Lagrangian solid models, which evaluate "
            << "d2dt2(D.oldTime()), are therefore not supported"
            << exit(FatalError);
    }

    const fvMesh& mesh = this->mesh();
    const Time& runTime = mesh.time();

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // The time at the start of the current time-step, to which the state
    // must be updated, and the field at that time. Before the first
    // time-step, for example when a solid model evaluates d2dt2 in its
    // constructor, this is the current time and field
    const scalar deltaT = runTime.deltaTValue();
    const bool firstStep = (runTime.timeIndex() == runTime.startTimeIndex());
    const scalar t0 = firstStep ? runTime.value() : runTime.value() - deltaT;

    // Request the old-time field on the first call too, as the Euler and
    // backward schemes do, so that it is stored from then on at the start of
    // each time-step, before the field is changed, for example by a boundary
    // condition or a fluid-solid interface update
    const fieldType& vfOld = vf.oldTime();
    const fieldType& vf0 = firstStep ? vf : vfOld;

    // Tolerance for comparing times
    const scalar tol = 1e-6*deltaT + 1e-12*mag(runTime.value());

    // The state time. What is written, and read on restart, is its lag
    // behind the time of writing, which keeps its precision when the times
    // are written with a low precision
    const word lagName("NewmarkLag(" + vfName + ')');
    bool newState = false;

    if (!mesh.foundObject<NewmarkLagField>(lagName))
    {
        const scalar startTime = runTime.startTime().value();

        IOobject lagIO
        (
            lagName,
            startTimeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

#ifdef OPENFOAM_NOT_EXTEND
        const bool present =
            lagIO.typeHeaderOk<uniformDimensionedScalarField>(true);
#else
        const bool present = lagIO.headerOk();
#endif

        scalar tState = t0;

        if (present)
        {
            const uniformDimensionedScalarField lag(lagIO);
            tState = startTime - lag.value();
        }
        else
        {
            newState = true;

            if (runTime.startTimeIndex() > 0)
            {
                WarningInFunction
                    << lagName << " was not found at the start time "
                    << startTimeName() << ": the " << type()
                    << " velocity and acceleration of " << vfName
                    << " are initialised from NewmarkV(" << vfName
                    << ") and NewmarkA(" << vfName << ") if present, and "
                    << "are zero otherwise" << endl;
            }

            if (t0 > startTime + tol)
            {
                WarningInFunction
                    << "The " << type() << " scheme is first evaluated for "
                    << vfName << " at time " << t0 << ", after the start "
                    << "time " << startTime << ": the initial velocity and "
                    << "acceleration are applied at time " << t0 << endl;
            }
        }

        regIOobject::store
        (
            new NewmarkLagField
            (
                IOobject
                (
                    lagName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                tState
            )
        );
    }

    scalar& stateTime =
        const_cast<NewmarkLagField&>
        (
            mesh.lookupObject<NewmarkLagField>(lagName)
        ).stateTime();

    fieldType& U =
        stateField
        (
            "NewmarkV(" + vfName + ')',
            dimensioned<Type>
            (
                "0", vf.dimensions()/dimTime, pTraits<Type>::zero
            )
        );

    fieldType& A =
        stateField
        (
            "NewmarkA(" + vfName + ')',
            dimensioned<Type>
            (
                "0", vf.dimensions()/dimTime/dimTime, pTraits<Type>::zero
            )
        );

    fieldType& D =
        stateField
        (
            "NewmarkD(" + vfName + ')',
            dimensioned<Type>("0", vf.dimensions(), pTraits<Type>::zero)
        );

    if (newState)
    {
        D == vf0;
    }

    scalar h = t0 - stateTime;

    // Before the first time-step, the state is not updated: when restarting,
    // it is updated in the first time-step, where deltaT0 is known
    if (h > tol && !firstStep)
    {
        // Complete the previous time-step, from the state time to t0, using
        // the converged field at t0
        const scalar deltaT0 = runTime.deltaT0Value();

        // Rounding error of a lag written with the write precision, and of
        // the difference of two times
        const scalar roundTol =
            0.5*Foam::pow(10.0, 1 - label(IOstream::defaultPrecision()))*h
          + 1e-14*mag(runTime.value());

        if (mag(h - deltaT0) <= roundTol)
        {
            // The interval is the previous time-step: use its exact value,
            // as the state time may have been read with a low precision
            h = deltaT0;
        }
        else
        {
            // The scheme was not evaluated in some time-steps, for example
            // because the solid is not solved before the fluid-solid
            // coupling starts, or the time-step was changed on restart. One
            // Newmark step is taken over the whole interval, which is exact
            // if the field did not change and its velocity and acceleration
            // were zero
            WarningInFunction
                << "The " << type() << " state of " << vfName << " is at "
                << "time " << stateTime << ", which is not one time-step ("
                << deltaT0 << ") before " << t0 << ": it is advanced to "
                << t0 << " in one step" << endl;
        }

        const dimensionedScalar hh("h", dimTime, h);

        const fieldType Anew
        (
            (vf0 - D - hh*U - (0.5 - beta_)*sqr(hh)*A)
           /(beta_*sqr(hh))
        );

        U == U + hh*((1.0 - gamma_)*A + gamma_*Anew);
        A == Anew;
        D == vf0;
        stateTime = t0;
    }
    else if (h < -tol)
    {
        FatalErrorInFunction
            << "The " << type() << " state of " << vfName << " is at time "
            << stateTime << ", which is after the start of the "
            << "current time-step at " << t0 << exit(FatalError);
    }
}


template<class Type>
void NewmarkBetaD2dt2Scheme<Type>::checkFvcScheme(const word& name) const
{
    // fvc::d2dt2 reads its scheme from ddtSchemes (issue #502), so, without
    // a ddtSchemes entry, an fvc::d2dt2 residual would silently use another
    // scheme than this fvm::d2dt2

#ifdef OPENFOAM_NOT_EXTEND
    ITstream& is = mesh().ddtScheme(name);
#else
    ITstream& is = mesh().schemesDict().ddtScheme(name);
#endif

    const word schemeName(is);

    bool same = (schemeName == typeName);

    if (same)
    {
        // Compare the coefficients
        const NewmarkBetaD2dt2Scheme<Type> other(mesh(), is);

        same =
            other.beta_ == beta_
         && other.gamma_ == gamma_
         && other.alphaM_ == alphaM_;
    }

    is.rewind();

    if (!same)
    {
        FatalErrorInFunction
            << "fvm::" << name << " uses the " << typeName
            << " scheme with beta = " << beta_ << ", gamma = " << gamma_
            << " and alphaM = " << alphaM_ << ", but fvc::" << name
            << ", which reads its scheme from ddtSchemes, would use a "
            << "different scheme or coefficients" << nl
            << "Add the same scheme to ddtSchemes, for example" << nl << nl
            << "    ddtSchemes" << nl
            << "    {" << nl
            << "        default         backward;" << nl
            << "        \"d2dt2\\(.*\\)\"   " << typeName << ";" << nl
            << "    }" << nl
            << exit(FatalError);
    }
}


template<class Type>
scalar NewmarkBetaD2dt2Scheme<Type>::diagCoeff() const
{
    return (1.0 - alphaM_)/(beta_*sqr(mesh().time().deltaTValue()));
}


template<class Type>
const GeometricField<Type, fvPatchField, volMesh>&
NewmarkBetaD2dt2Scheme<Type>::state
(
    const word& prefix,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return mesh().template
        lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            prefix + '(' + vf.name() + ')'
        );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::explicitPart
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    updateState(vf);

    return explicitPart
    (
        vf.name(),
        state("NewmarkV", vf),
        state("NewmarkA", vf),
        state("NewmarkD", vf)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::explicitPart
(
    const word& vfName,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvPatchField, volMesh>& A,
    const GeometricField<Type, fvPatchField, volMesh>& D
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const dimensionedScalar deltaT = mesh().time().deltaT();
    const dimensionedScalar c
    (
        "c", dimless/dimTime/dimTime, diagCoeff()
    );

    return tmp<fieldType>
    (
        new fieldType
        (
            IOobject
            (
                "d2dt2Explicit(" + vfName + ')',
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            c*(D + deltaT*U + (0.5 - beta_)*sqr(deltaT)*A) - alphaM_*A
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::acceleration
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    updateState(vf);

    const fieldType& A = state("NewmarkA", vf);

    const Time& runTime = mesh().time();

    if (runTime.timeIndex() == runTime.startTimeIndex())
    {
        // Before the first time-step, return the initial acceleration
        return tmp<fieldType>(new fieldType(A));
    }

    const dimensionedScalar c
    (
        "c", dimless/dimTime/dimTime, diagCoeff()
    );

    return
        c*vf
      - explicitPart
        (
            vf.name(), state("NewmarkV", vf), A, state("NewmarkD", vf)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
NewmarkBetaD2dt2Scheme<Type>::NewmarkBetaD2dt2Scheme
(
    const fvMesh& mesh,
    Istream& is
)
:
    d2dt2Scheme<Type>(mesh, is),
    beta_(0.25),
    gamma_(0.5),
    alphaM_(0.0)
{
    // Optional coefficients: beta, gamma and alphaM
    if (!is.eof())
    {
        beta_ = readScalar(is);
    }

    if (!is.eof())
    {
        gamma_ = readScalar(is);
    }

    if (!is.eof())
    {
        alphaM_ = readScalar(is);
    }

    if (beta_ <= 0 || gamma_ <= 0 || alphaM_ >= 1)
    {
        FatalIOErrorInFunction(is)
            << "The " << type() << " coefficients must satisfy beta > 0, "
            << "gamma > 0 and alphaM < 1, but they are beta = " << beta_
            << ", gamma = " << gamma_ << " and alphaM = " << alphaM_
            << exit(FatalIOError);
    }

    // Warn about coefficients outside the documented ranges, where the
    // scheme is not unconditionally stable or adds negative damping. The
    // scheme is constructed for every fvm::d2dt2 and fvc::d2dt2 call, so each
    // warning is given once
    static bool warnedGamma = false;
    static bool warnedBeta = false;
    static bool warnedAlphaM = false;

    if (gamma_ < 0.5 && !warnedGamma)
    {
        warnedGamma = true;

        WarningInFunction
            << "gamma = " << gamma_ << " < 1/2: the " << type()
            << " scheme amplifies the motion (negative numerical damping)"
            << endl;
    }

    if (beta_ < 0.5*gamma_ && !warnedBeta)
    {
        warnedBeta = true;

        WarningInFunction
            << "beta = " << beta_ << " < gamma/2 = " << 0.5*gamma_
            << ": the " << type() << " scheme is only conditionally stable"
            << endl;
    }

    if ((alphaM_ < -1.0/3.0 || alphaM_ > 0) && !warnedAlphaM)
    {
        warnedAlphaM = true;

        WarningInFunction
            << "alphaM = " << alphaM_ << " is outside the Bossak range "
            << "-1/3 <= alphaM <= 0" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::physicalAcceleration
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    tmp<fieldType> ta = fvc::d2dt2(vf);

    // fvc::d2dt2 reads its scheme from ddtSchemes (issue #502)
    const fvMesh& mesh = vf.mesh();
    const word name("d2dt2(" + vf.name() + ')');

#ifdef OPENFOAM_NOT_EXTEND
    ITstream& is = mesh.ddtScheme(name);
#else
    ITstream& is = mesh.schemesDict().ddtScheme(name);
#endif

    const word schemeName(is);

    scalar alphaM = 0;

    if (schemeName == typeName)
    {
        // The optional coefficients: beta, gamma and alphaM
        for (label coeffI = 0; coeffI < 3 && !is.eof(); ++coeffI)
        {
            const scalar coeff = readScalar(is);

            if (coeffI == 2)
            {
                alphaM = coeff;
            }
        }
    }

    is.rewind();

    if (alphaM == 0)
    {
        return ta;
    }

    // d2dt2 is the Bossak-weighted (1 - alphaM) a^{n+1} + alphaM a^n, and
    // the stored NewmarkA is a^n once the state has been updated in this
    // time-step, as fvc::d2dt2 did above; before the first time-step, both
    // are the initial acceleration
    const fieldType& A =
        mesh.template lookupObject<fieldType>("NewmarkA(" + vf.name() + ')');

    return (ta - alphaM*A)/(1.0 - alphaM);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tacc = acceleration(vf);

    tmpRef(tacc).rename("d2dt2(" + vf.name() + ')');

    return tacc;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tacc =
        rho*acceleration(vf);

    tmpRef(tacc).rename("d2dt2(" + rho.name() + ',' + vf.name() + ')');

    return tacc;
}


template<class Type>
tmp<fvMatrix<Type> >
NewmarkBetaD2dt2Scheme<Type>::fvmD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tmpRef(tfvm);

    checkFvcScheme("d2dt2(" + vf.name() + ')');

    const scalarField& V = mesh().V();
    const tmp<GeometricField<Type, fvPatchField, volMesh> > texplicit =
        explicitPart(vf);

    fvm.diag() = diagCoeff()*V;
    fvm.source() = V*primitiveField(texplicit());

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
NewmarkBetaD2dt2Scheme<Type>::fvmD2dt2
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
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tmpRef(tfvm);

    checkFvcScheme("d2dt2(" + vf.name() + ')');

    const scalarField& V = mesh().V();
    const tmp<GeometricField<Type, fvPatchField, volMesh> > texplicit =
        explicitPart(vf);

    fvm.diag() = diagCoeff()*rho.value()*V;
    fvm.source() = rho.value()*V*primitiveField(texplicit());

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
NewmarkBetaD2dt2Scheme<Type>::fvmD2dt2
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
            vf.dimensions()*rho.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tmpRef(tfvm);

    checkFvcScheme("d2dt2(" + rho.name() + ',' + vf.name() + ')');

    const scalarField& V = mesh().V();
    const tmp<GeometricField<Type, fvPatchField, volMesh> > texplicit =
        explicitPart(vf);

    // The current density multiplies the acceleration
    const scalarField& rhoI = primitiveField(rho);
    fvm.diag() = diagCoeff()*rhoI*V;
    fvm.source() = rhoI*V*primitiveField(texplicit());

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
