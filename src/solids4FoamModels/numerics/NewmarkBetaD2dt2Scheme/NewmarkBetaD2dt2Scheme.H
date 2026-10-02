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

Class
    Foam::fv::NewmarkBetaD2dt2Scheme

Description
    Newmark-beta d2dt2 scheme, with optional Bossak-alpha damping.

    For the time-step from t^n to t^{n+1} = t^n + deltaT, the acceleration is

        a^{n+1} = (u^{n+1} - ubar)/(beta*deltaT^2)
        ubar = u^n + deltaT*v^n + (0.5 - beta)*deltaT^2*a^n

    and the scheme returns the Bossak-weighted acceleration

        d2dt2(u) = (1 - alphaM)*a^{n+1} + alphaM*a^n

    The velocity and acceleration are then updated as

        v^{n+1} = v^n + deltaT*((1 - gamma)*a^n + gamma*a^{n+1})

    Usage in system/fvSchemes:

        d2dt2Schemes
        {
            default         NewmarkBeta;  // beta = 1/4, gamma = 1/2
        }

        ddtSchemes
        {
            default         backward;
            "d2dt2\(.*\)"   NewmarkBeta;
        }

    The optional coefficients follow the name: "NewmarkBeta beta gamma" or
    "NewmarkBeta beta gamma alphaM", with the defaults beta = 0.25,
    gamma = 0.5 and alphaM = 0.

    - beta = 1/4, gamma = 1/2, alphaM = 0 is the average acceleration
      (trapezoidal) rule: second-order, unconditionally stable and without
      numerical damping.
    - gamma > 1/2 with beta = (gamma + 1/2)^2/4 adds damping, but the scheme
      is then first-order.
    - The Bossak-alpha method uses -1/3 <= alphaM <= 0,
      gamma = 1/2 - alphaM and beta = (1 - alphaM)^2/4. It is second-order
      and damps high frequencies, with the spectral radius at infinity
      rho_inf = (1 + alphaM)/(1 - alphaM).

    The ddtSchemes entry is needed because fvc::d2dt2, which is used for
    example in the PETSc SNES residual, reads its scheme from ddtSchemes
    (issue #502). Without it, the residual would silently use another scheme,
    so fvm::d2dt2 raises a fatal error if the ddtSchemes entry for d2dt2(D)
    is not NewmarkBeta. There is no NewmarkBeta ddt scheme, so it must not be
    the ddtSchemes default. The velocity from fvc::ddt keeps the ddtSchemes
    default scheme, and it is not the Newmark velocity.

    State: the velocity, the acceleration and the field at the last completed
    time-step are stored as the registered fields "NewmarkV(D)",
    "NewmarkA(D)" and "NewmarkD(D)", where D is the field name, and the lag
    of that time-step behind the time of writing as "NewmarkLag(D)". They are
    written with the other fields, and they are read if present, so a
    restart continues the time integration; the lag, rather than the time, is
    written so that it keeps its precision. They are updated once per
    time-step, on the first call in a new time-step, from the old-time field,
    which is then the converged field of the completed time-step. The
    coupling and nonlinear iterations within a time-step therefore do not
    change the state. If the scheme is not evaluated for some time-steps, for
    example before a fluid-solid coupling start time, the state is advanced
    over that interval in one step, with a warning; this is exact if the
    field did not change and its velocity and acceleration were zero.

    Initial conditions: at the start of a run, the velocity and the
    acceleration are read from "NewmarkV(D)" and "NewmarkA(D)" in the start
    time directory if they are present, and are zero otherwise. A zero
    initial acceleration assumes that the body starts in equilibrium with the
    initial loads, for example at rest and unloaded. Otherwise, provide the
    initial acceleration, for example from the equilibrium equation, in
    "NewmarkA(D)". Before the first time-step, fvc::d2dt2 returns the stored
    acceleration: the initial acceleration, or, on restart, that of the time
    before the last completed time-step.

    With alphaM != 0, fvc::d2dt2 returns the Bossak-weighted acceleration
    (1 - alphaM) a^{n+1} + alphaM a^n, which the momentum residual needs.
    Where the physical acceleration a^{n+1} is meant, such as the interface
    acceleration passed to the fluid (solidModel::faceZoneAcceleration) and
    the acceleration of the total Lagrangian solid model, it is recovered
    with physicalAcceleration().

    The time-step may vary. The scheme is not implemented for a moving mesh,
    or for the updated Lagrangian solid models, which evaluate
    d2dt2(D.oldTime()): a fatal error is raised for a field whose name ends
    in "_0". The vertex-centred solid models do not use this scheme: they
    have their own NewmarkBeta implementation, which requires beta and gamma
    in d2dt2Schemes and NewmarkBeta in ddtSchemes.

SourceFiles
    NewmarkBetaD2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef NewmarkBetaD2dt2Scheme_H
#define NewmarkBetaD2dt2Scheme_H

#include "d2dt2Scheme.H"
#include "uniformDimensionedFields.H"
#ifdef FOAMEXTEND
    #include "foamTime.H"
#else
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class NewmarkLagField Declaration
\*---------------------------------------------------------------------------*/

//- The time of the NewmarkBeta state, which is written as its lag behind
//  the time of writing
class NewmarkLagField
:
    public uniformDimensionedScalarField
{
    // Private Data

        //- The state time
        scalar stateTime_;


public:

    // Constructors

        //- Construct from IOobject and the state time
        NewmarkLagField(const IOobject& io, const scalar stateTime)
        :
            uniformDimensionedScalarField
            (
                io,
                dimensionedScalar("lag", dimTime, 0)
            ),
            stateTime_(stateTime)
        {}


    // Member Functions

        //- Return the state time
        scalar& stateTime()
        {
            return stateTime_;
        }

        //- Write the lag of the state time behind the current time
        virtual bool writeData(Ostream& os) const
        {
            const_cast<NewmarkLagField&>(*this).value() =
                time().value() - stateTime_;

            return uniformDimensionedScalarField::writeData(os);
        }
};


/*---------------------------------------------------------------------------*\
                   Class NewmarkBetaD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NewmarkBetaD2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
    // Private Data

        //- Newmark beta coefficient
        scalar beta_;

        //- Newmark gamma coefficient
        scalar gamma_;

        //- Bossak alpha coefficient
        scalar alphaM_;


    // Private Member Functions

        //- Return the name of the start time directory
        word startTimeName() const;

        //- Return the named state field, creating it if needed: it is read
        //  if present, and is otherwise initialised to the given value
        GeometricField<Type, fvPatchField, volMesh>& stateField
        (
            const word& name,
            const dimensioned<Type>& value
        ) const;

        //- Update the state to the start of the current time-step, if not
        //  already done
        void updateState
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Check that fvc::d2dt2, which reads the named scheme from
        //  ddtSchemes, uses this scheme with the same coefficients
        void checkFvcScheme(const word& name) const;

        //- Return the coefficient of u^{n+1} in d2dt2(u)
        scalar diagCoeff() const;

        //- Return the state field prefix(u), e.g. NewmarkA(u), which must
        //  exist
        const GeometricField<Type, fvPatchField, volMesh>& state
        (
            const word& prefix,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return the explicit part of d2dt2(u), so that
        //  d2dt2(u) = diagCoeff()*u^{n+1} - explicitPart(u), after updating
        //  the state
        tmp<GeometricField<Type, fvPatchField, volMesh> > explicitPart
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return the explicit part of d2dt2(u) from the given state: the
        //  velocity U, the acceleration A and the field D
        tmp<GeometricField<Type, fvPatchField, volMesh> > explicitPart
        (
            const word& vfName,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvPatchField, volMesh>& A,
            const GeometricField<Type, fvPatchField, volMesh>& D
        ) const;

        //- Return d2dt2(u), or the initial acceleration before the first
        //  time-step
        tmp<GeometricField<Type, fvPatchField, volMesh> > acceleration
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Disallow default bitwise copy construct
        NewmarkBetaD2dt2Scheme(const NewmarkBetaD2dt2Scheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const NewmarkBetaD2dt2Scheme&) = delete;


public:

    //- Runtime type information
    TypeName("NewmarkBeta");


    // Constructors

        //- Construct from mesh and Istream
        NewmarkBetaD2dt2Scheme(const fvMesh& mesh, Istream& is);


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

        //- Return the physical acceleration of vf, a^{n+1}. With the
        //  Bossak-alpha method (alphaM != 0), fvc::d2dt2 returns the weighted
        //  (1 - alphaM) a^{n+1} + alphaM a^n, which the momentum residual
        //  needs; this recovers a^{n+1} from it. For any other d2dt2 scheme,
        //  it returns fvc::d2dt2(vf)
        static tmp<GeometricField<Type, fvPatchField, volMesh> >
        physicalAcceleration
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NewmarkBetaD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
