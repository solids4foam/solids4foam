/*---------------------------------------------------------------------------*\
  XenoSim - Cardiovascular CFD & Multiphysics
  Component: BDFDdtScheme

  Author: David Kelly
  Copyright (C) 2026 David Kelly

Class
    Foam::fv::BDFDdtScheme

Description
    Constant-step backward differentiation formula (BDF) ddt of a fixed
    maximum order: the shared implementation behind the BDF3 and BDF4
    schemes, which differ only in that maximum.

    \verbatim
        ddt(phi)^(n+1) = (1/dt) sum_{k=0..q} a_k phi^(n+1-k)
    \endverbatim

    with the coefficients hard-coded for each order q:

    \table
        q | a_0   | a_1 | a_2 | a_3  | a_4
        1 | 1     | -1  |     |      |
        2 | 3/2   | -2  | 1/2 |      |
        3 | 11/6  | -3  | 3/2 | -1/3 |
        4 | 25/12 | -4  | 3   | -4/3 | 1/4
    \endtable

    Start-up. A q-th order formula needs q old time levels. Each level's
    time index records the step that filled it; a level created on demand is
    a copy of the one above and shares its index. The order used on a step is
    therefore one more than the number of consecutive distinct pairs of old
    levels, up to the maximum, so a fresh start ramps Euler, BDF2, BDF3, BDF4
    over its first steps, exactly as OpenFOAM's backward scheme falls back to
    Euler on its first step. Those start-up steps carry lower-order local
    errors, so to see the full order in a convergence study either supply the
    history - write the old-time fields (U_0, U_0_0, ...) into the start time,
    which OpenFOAM reads back with distinct time indices - or start from a
    solution whose low derivatives vanish at t = 0. A restart finds the
    old-time fields the previous run wrote and resumes at full order.

    Every old level up to the maximum is created on the first step, whatever
    order that step uses, so that each is filled as the run advances. Creating
    a level only when the order first needs it would fill it with a copy of
    the level above instead of the genuine older value.

    Time step. The coefficients are for a uniform step. Time records only the
    current and previous step, so a change of step is detected on the step it
    happens and is fatal, rather than silently applying uniform weights to a
    non-uniform stencil. A restart with a changed deltaT is not detectable and
    must be avoided.

    Moving meshes are not supported: fvMesh keeps only V0 and V00, fewer old
    volumes than a third- or fourth-order formula needs, so construction on a
    moving mesh is fatal.

    The flux corrections used by the pressure-velocity coupling (ddtCorr) are
    the same q-th order combination applied to the old fluxes and the old
    interpolated velocities, the direct generalisation of backward's.

SourceFiles
    BDFDdtScheme.C
    BDFDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef BDFDdtScheme_H
#define BDFDdtScheme_H

#include "ddtScheme.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class BDFDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BDFDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Data

        //- Order used once enough old time levels are available
        const label maxOrder_;


    // Private Member Functions

        //- Old time level k of vf; k = 0 is vf itself. Levels are created on
        //  demand, and accessing level 1 stores the old times when the time
        //  index has moved on.
        template<class GeoField>
        static const GeoField& level(const GeoField& vf, const label k);

        //- Create every old level of vf up to maxOrder_
        template<class GeoField>
        void touch(const GeoField& vf) const;

        //- Order usable on this step for vf: one more than the number of
        //  consecutive old levels holding distinct time indices
        template<class GeoField>
        label order(const GeoField& vf) const;

        //- Difference weights a_0..a_q for order q. Fails if q > 1 and the
        //  time step has changed.
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

        //- Fail on a moving mesh
        void checkStaticMesh() const;

        //- IOobject for a returned ddt field
        IOobject ddtIOobject(const word& name) const;

        //- No copy construct
        BDFDdtScheme(const BDFDdtScheme&) = delete;

        //- No copy assignment
        void operator=(const BDFDdtScheme&) = delete;


protected:

    // Constructors

        //- Construct from mesh with the given maximum order
        BDFDdtScheme(const fvMesh& mesh, const label maxOrder);

        //- Construct from mesh and Istream with the given maximum order. An
        //  optional ddtPhiCoeff may follow the scheme name, as for backward.
        BDFDdtScheme(const fvMesh& mesh, Istream& is, const label maxOrder);


public:

    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        //- Maximum order
        label maxOrder() const noexcept
        {
            return maxOrder_;
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
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtPhiCorr
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
    #include "BDFDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
