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
    Foam::fv::bdfDdtScheme

Description
    Implicit backward differentiation formula (BDF) ddt of order 1 to 6. BDF 1
    and BDF 2 allow a variable time-step; BDF 3 to BDF 6 require a constant
    one.

    Usage in system/fvSchemes, where the order is 1 to 6:

        ddtSchemes
        {
            default BDF 3;
        }

    BDF 1 is the Euler scheme and BDF 2 is the backward scheme, for a constant
    or a variable time-step. This is a ddt scheme, for fvm::ddt and
    fvc::ddt; for the second time derivative, see bdfD2dt2Scheme.

    \verbatim
        ddt(phi)^(n+1) = (1/dt) sum_{k=0..q} a_k phi^(n+1-k)
    \endverbatim

    with the coefficients hard-coded for each order q:

    \table
        q | a_0    | a_1 | a_2  | a_3   | a_4  | a_5  | a_6
        1 | 1      | -1  |      |       |      |      |
        2 | 3/2    | -2  | 1/2  |       |      |      |
        3 | 11/6   | -3  | 3/2  | -1/3  |      |      |
        4 | 25/12  | -4  | 3    | -4/3  | 1/4  |      |
        5 | 137/60 | -5  | 5    | -10/3 | 5/4  | -1/5 |
        6 | 49/20  | -6  | 15/2 | -20/3 | 15/4 | -6/5 | 1/6
    \endtable

    Warning: BDF 3 to BDF 6 are not A-stable, only A(alpha)-stable: their
    stability regions contain the sector |arg(-z)| < alpha but leave out a
    region around the imaginary axis, which widens with the order:

    \table
        q     | 3     | 4     | 5     | 6
        alpha | 86.0° | 73.4° | 51.8° | 17.8°
    \endtable

    BDF 7 and above are zero-unstable, hence the cap at 6. BDF 3 and 4 are
    standard for diffusion-dominated, stiff problems. In strongly
    convection-dominated flow at large Courant numbers they can be less robust
    than backward, and BDF 5 and 6 are not recommended for it. A warning is
    printed when an order of 3 or above is selected.

    Start-up. A q-th order formula needs q old time levels. Each level's time
    index records the step that filled it; a level created on demand is a copy
    of the one above and shares its index. The order used on a step is
    therefore one more than the number of consecutive distinct pairs of old
    levels, up to the maximum, so a fresh start ramps Euler, BDF2, BDF3, ...
    up to the selected order over its first steps, exactly as OpenFOAM's
    backward scheme falls back to Euler on its first step. Those start-up
    steps carry lower-order local errors, so to see the full order in a
    convergence study either supply the history - write the old-time fields
    (U_0, U_0_0, ...) into the start time, which OpenFOAM reads back with
    distinct time indices - or start from a solution whose low derivatives
    vanish at t = 0. A restart finds the old-time fields the previous run
    wrote and resumes at full order.

    Fields without stored history. A field may have fewer stored old levels
    than U, e.g. a rho derived rather than read on restart. Its copied levels
    must not be weighted as genuine older values:

    - rho and alpha, in ddt(rho, U) and ddt(alpha, rho, U) and in the flux
      corrections: the order is the lowest available across all of them, so
      it drops and ramps up again as their levels fill.
    - The old fluxes phi or Uf of the flux corrections: the order stays that
      of U, matching the momentum equation, and each flux level older than
      the stored history is reconstructed from the same level of U (or
      rho*U), so that it contributes nothing to the correction.

    Either way the first step a field exists on uses a copy of its current
    value as its previous one, as no older value is known.

    Every old level up to the maximum is created on the first step, whatever
    order that step uses, so that each is filled as the run advances. Creating
    a level only when the order first needs it would fill it with a copy of
    the level above instead of the genuine older value.

    Cost. Each field the scheme is applied to, and the fluxes of the flux
    corrections, therefore keep q old levels in memory from the first step,
    q + 1 copies of the field in all. At each write time q - 1 of them are
    written as well, which is all a restart needs: BDF 3 writes U_0 and
    U_0_0, and BDF 6 writes U_0 to U_0_0_0_0_0.

    Time step. BDF 1 and BDF 2 support a variable time step, including
    adjustTimeStep yes. BDF 2 then uses the variable-step coefficients

    \verbatim
        a_0 = (1 + 2 w)/(1 + w),  a_1 = -(1 + w),  a_2 = w^2/(1 + w)
    \endverbatim

    with w = deltaT/deltaT0, which are exactly those of OpenFOAM's backward
    scheme, so BDF 2 is a drop-in replacement for backward on a static mesh.

    BDF 3 to BDF 6 require a constant time step, so adjustTimeStep yes cannot
    be used with them. Their variable-step coefficients would need the older
    step sizes as well, but Time records only the current and previous step.
    A change of step is therefore detected on the step it happens and is
    fatal, rather than silently applying uniform weights to a non-uniform
    stencil.

    Restart with a changed deltaT. With adjustTimeStep no, Time takes deltaT
    from controlDict on restart and sets deltaT0 to it on the first step, so
    the step the old levels were stored with is lost. On that first step the
    scheme therefore reads it from <startTime>/uniform/time instead: BDF 2
    uses it in its variable-step weights, and BDF 3 to BDF 6 stop with a
    fatal error if it differs from the new deltaT. A change is detected to a
    relative tolerance of 1e-5, as uniform/time holds deltaT to writePrecision
    digits only.

    Moving meshes are not supported: fvMesh keeps only V0 and V00, fewer old
    volumes than a third- or higher-order formula needs, so construction on a
    moving mesh is fatal.

    The flux corrections used by the pressure-velocity coupling (ddtCorr) are
    the same q-th order combination applied to the old fluxes and the old
    interpolated velocities, the direct generalisation of backward's.

Author
    David Kelly, UCD.

SourceFiles
    bdfDdtScheme.C
    bdfDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef bdfDdtScheme_H
#define bdfDdtScheme_H

#include "ddtScheme.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

//- Whether the order 3 and above warning has been printed. It is shared by
//  every instantiation, so the warning is printed once per run, not once per
//  field type.
extern bool bdfDdtSchemeWarned;


/*---------------------------------------------------------------------------*\
                        Class bdfDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bdfDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Data

        //- Order used once enough old time levels are available
        const label order_;

        //- Relative tolerance on a change of time step, as in
        //  bdfD2dt2Scheme. It allows for a deltaT read back from
        //  uniform/time, which has only writePrecision (default 6) digits.
        static constexpr scalar deltaTTol_ = 1e-5;


    // Private Member Functions

        //- Old time level k of vf; k = 0 is vf itself. Levels are created on
        //  demand, and accessing level 1 stores the old times when the time
        //  index has moved on.
        template<class GeoField>
        static const GeoField& level(const GeoField& vf, const label k);

        //- Create every old level of vf up to order_
        template<class GeoField>
        void touch(const GeoField& vf) const;

        //- Order usable on this step for vf: one more than the number of
        //  consecutive old levels holding distinct time indices
        template<class GeoField>
        label order(const GeoField& vf) const;

        //- Order usable on this step for a formula involving two fields:
        //  the lower of their orders
        template<class GeoField1, class GeoField2>
        label order(const GeoField1& f1, const GeoField2& f2) const;

        //- Order usable on this step for a formula involving three fields:
        //  the lowest of their orders
        template<class GeoField1, class GeoField2, class GeoField3>
        label order
        (
            const GeoField1& f1,
            const GeoField2& f2,
            const GeoField3& f3
        ) const;

        //- The previous time step. On the first step after a restart this
        //  is the deltaT stored in <startTime>/uniform/time, if it differs
        //  from Time's deltaT0, which is then the new deltaT of controlDict.
        scalar deltaT0() const;

        //- Difference weights a_0..a_q for order q. BDF 2 uses the
        //  variable-step weights; otherwise fails if q > 1 and the time step
        //  has changed.
        scalarList weights(const label q) const;

        //- The old-time part of the formula, sum_{k=1..q} (-a_k) vf^(n+1-k)
        template<class GeoField>
        tmp<GeoField> history(const GeoField& vf, const scalarList& a) const;

        //- As history(vf), for the product rho*vf
        tmp<GeometricField<Type, fvPatchField, volMesh>> history
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const scalarList& a
        ) const;

        //- As history(vf), for the product alpha*rho*vf
        tmp<GeometricField<Type, fvPatchField, volMesh>> history
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const scalarList& a
        ) const;

        //- As history(flux), for the old fluxes phi or Uf of the flux
        //  corrections, but each level older than the stored flux history
        //  is replaced by reconstruct(k), the flux interpolated from level
        //  k of the cell field, which contributes nothing to the correction
        template<class GeoField, class Reconstruct>
        tmp<GeoField> fluxHistory
        (
            const GeoField& flux,
            const scalarList& a,
            const Reconstruct& reconstruct
        ) const;

        //- Fail on a moving mesh
        void checkStaticMesh() const;

        //- IOobject for a returned ddt field
        IOobject ddtIOobject(const word& name) const;

        //- Read and validate the BDF order from an Istream, raising a
        //  fatal error if it is missing or not 1 to 6
        static label readOrder(Istream& is);

        //- No copy construct
        bdfDdtScheme(const bdfDdtScheme&) = delete;

        //- No copy assignment
        void operator=(const bdfDdtScheme&) = delete;


public:

    //- Runtime type information
    TypeName("BDF");


    // Constructors

        //- Construct from mesh and Istream. The order follows the scheme
        //  name, e.g. "BDF 3"; an optional ddtPhiCoeff may follow the order,
        //  as for backward.
        bdfDdtScheme(const fvMesh& mesh, Istream& is);


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// The flux corrections have no meaning for a scalar "velocity"; these are
// specialised to fail, as they are for every OpenFOAM ddt scheme.

template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "bdfDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
