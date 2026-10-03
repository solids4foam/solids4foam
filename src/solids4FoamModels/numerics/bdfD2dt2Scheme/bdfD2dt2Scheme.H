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
    Foam::fv::bdfD2dt2Scheme

Description
    Implicit backward differentiation formula (BDF) d2dt2 of order 1 to 6,
    for a constant time-step.

    The BDFk first derivative is applied twice, in the same way as the
    backward d2dt2 scheme applies the backward ddt scheme twice:

        v^n = (1/deltaT) sum_{i=0}^{k} a_i u^{n-i}
        a^n = (1/deltaT) sum_{j=0}^{k} a_j v^{n-j}

    where a_i are the BDFk coefficients. This is BDFk applied to the
    first-order system du/dt = v, dv/dt = a, so it inherits the stability of
    BDFk. It uses the current and 2k previous time levels. The velocity
    implied by the scheme is recomputed from the stored displacements, so no
    velocity field is stored.

    Usage in system/fvSchemes, where the order is 1 to 6:

        d2dt2Schemes
        {
            default BDF 3;
        }

        ddtSchemes
        {
            default         backward;
            "d2dt2\(.*\)"   BDF 3;
        }

    The ddtSchemes entry is needed because fvc::d2dt2 in OpenFOAM and
    foam-extend reads its scheme from ddtSchemes, while fvm::d2dt2 reads it
    from d2dt2Schemes. Without it, fvc::d2dt2, which is used for example in
    the PETSc SNES residual, falls back to the ddtSchemes default. The
    velocity from fvc::ddt keeps the ddtSchemes default scheme. The
    vertex-centred solid models select their own d2dt2 schemes and do not use
    this scheme.

    With a constant time-step and a constant density, "BDF 1" is identical
    to Euler, and, once 2k old-time levels are available, "BDF 2" is
    identical to the solids4foam backward d2dt2 scheme.

    Start-up: when fewer than 2k genuine old-time levels are available, the
    order is reduced to half of the number of available levels (at least
    one), as the backward scheme falls back to Euler. A missing level is
    detected as an old-time field with the same time index as its parent, as
    in the backward scheme; on OpenFOAM.org and foam-extend this needs the
    GeometricField.C fix in optionalFixes. The old-time levels are written
    with the field (with "restart yes" in the solid model coefficients), so a
    restart continues at full order. A restart within the first 2k time-steps
    can read copies of the initial field as genuine levels, so its order can
    rise earlier than in a continuous run, and its results can differ.

    The time-step must be constant: a fatal error is raised if the time-step
    changes while two or more genuine old-time levels are available. The
    scheme is not implemented for a moving mesh.

    Note: BDF3 to BDF6 are not A-stable. For an undamped oscillator with
    angular frequency omega, the BDF3 amplitude grows for omega*deltaT below
    1.94 and the BDF4 amplitude grows for omega*deltaT below 4.71. At 20
    steps per period the growth is about 4.5% per period for BDF3 and 0.6%
    per period for BDF4. BDF5 and BDF6 slightly damp well-resolved modes,
    but a spurious root grows, by up to 38% and 58% per time-step, for
    omega*deltaT between 0.71 and 9.39 (BDF5) and between 0.84 and 17.6
    (BDF6), that is for modes with fewer than about 9 and 7.5 steps per
    period. A mesh normally has such modes, so BDF5 and BDF6 are unstable
    unless the time-step resolves the highest mesh frequency. BDF3 to BDF6
    are therefore not recommended for undamped structural dynamics;
    "BDF 1" and "BDF 2" are dissipative.

SourceFiles
    bdfD2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef bdfD2dt2Scheme_H
#define bdfD2dt2Scheme_H

#include "d2dt2Scheme.H"
#include "scalarList.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class bdfD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bdfD2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
    // Private Data

        //- Order of the scheme
        const label order_;


    // Private Member Functions

        //- Return the BDF first derivative coefficients of the given order
        static scalarList coeffs(const label order);

        //- Return the n-th old-time level of the field
        template<class GeoField>
        static const GeoField& oldTimeLevel
        (
            const GeoField& vf,
            const label n
        );

        //- Return the order to use this time-step: the order is reduced
        //  during start-up, and a fatal error is raised if the time-step
        //  has changed
        label order
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return sum_{i=i0}^{p} a_i u^{n-j-i}, where a_i are the BDF
        //  coefficients of order p and j is the old-time level
        tmp<Field<Type> > levelSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const label j,
            const label i0,
            const label p
        ) const;

        //- Return the BDF first derivative of order p at old-time level j
        tmp<GeometricField<Type, fvPatchField, volMesh> > levelDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const label j,
            const label p
        ) const;

        //- Check that the mesh is not moving
        void checkMesh() const;

        //- Disallow default bitwise copy construct
        bdfD2dt2Scheme(const bdfD2dt2Scheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const bdfD2dt2Scheme&) = delete;


public:

    //- Runtime type information
    TypeName("BDF");


    // Constructors

        //- Construct from mesh and Istream
        bdfD2dt2Scheme(const fvMesh& mesh, Istream& is);


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

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
#   include "bdfD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
