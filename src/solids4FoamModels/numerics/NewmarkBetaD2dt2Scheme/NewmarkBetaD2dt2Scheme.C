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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
        mesh.time().timeName(mesh.time().startTime().value()),
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
    const fieldType& vf0 = firstStep ? vf : vf.oldTime();

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
            runTime.timeName(startTime),
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
                    << runTime.timeName(startTime) << ": the " << type()
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
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::explicitPart
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    updateState(vf);

    const word& vfName = vf.name();
    const objectRegistry& db = mesh();

    const fieldType& U =
        db.lookupObject<fieldType>("NewmarkV(" + vfName + ')');
    const fieldType& A =
        db.lookupObject<fieldType>("NewmarkA(" + vfName + ')');
    const fieldType& D =
        db.lookupObject<fieldType>("NewmarkD(" + vfName + ')');

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

    const tmp<fieldType> texplicit = explicitPart(vf);

    const Time& runTime = mesh().time();

    if (runTime.timeIndex() == runTime.startTimeIndex())
    {
        // Before the first time-step, return the initial acceleration
        const objectRegistry& db = mesh();

        return tmp<fieldType>
        (
            new fieldType
            (
                db.lookupObject<fieldType>("NewmarkA(" + vf.name() + ')')
            )
        );
    }

    const dimensionedScalar c
    (
        "c", dimless/dimTime/dimTime, diagCoeff()
    );

    return c*vf - texplicit;
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh> >
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
            acceleration(vf)
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
NewmarkBetaD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh> >
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
            rho*acceleration(vf)
        )
    );
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

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    checkFvcScheme("d2dt2(" + vf.name() + ')');

    const scalarField& V = mesh().V();
    const tmp<GeometricField<Type, fvPatchField, volMesh> > texplicit =
        explicitPart(vf);

    fvm.diag() = diagCoeff()*V;
#ifdef FOAMEXTEND
    fvm.source() = V*texplicit().internalField();
#else
    fvm.source() = V*texplicit().primitiveField();
#endif

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

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    checkFvcScheme("d2dt2(" + vf.name() + ')');

    const scalarField& V = mesh().V();
    const tmp<GeometricField<Type, fvPatchField, volMesh> > texplicit =
        explicitPart(vf);

    fvm.diag() = diagCoeff()*rho.value()*V;
#ifdef FOAMEXTEND
    fvm.source() = rho.value()*V*texplicit().internalField();
#else
    fvm.source() = rho.value()*V*texplicit().primitiveField();
#endif

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

#ifdef FOAMEXTEND
    fvMatrix<Type>& fvm = tfvm();
#else
    fvMatrix<Type>& fvm = tfvm.ref();
#endif

    checkFvcScheme("d2dt2(" + rho.name() + ',' + vf.name() + ')');

    const scalarField& V = mesh().V();
    const tmp<GeometricField<Type, fvPatchField, volMesh> > texplicit =
        explicitPart(vf);

    // The current density multiplies the acceleration
#ifdef FOAMEXTEND
    const scalarField& rhoI = rho.internalField();
    fvm.diag() = diagCoeff()*rhoI*V;
    fvm.source() = rhoI*V*texplicit().internalField();
#else
    const scalarField& rhoI = rho.primitiveField();
    fvm.diag() = diagCoeff()*rhoI*V;
    fvm.source() = rhoI*V*texplicit().primitiveField();
#endif

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
