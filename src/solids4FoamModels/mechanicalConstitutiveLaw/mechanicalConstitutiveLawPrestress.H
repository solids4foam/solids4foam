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
    Foam::mechanicalConstitutiveLawPrestress

Description
    The prestress of a mechanicalConstitutiveLawManager, made when
    "prestress yes;" is given in mechanicalProperties: the prestress
    deformation gradient Fpre, read from the start time, which maps a
    stress-free configuration onto the mesh the solver is solving on, with
    its inverse and determinant, and the same interpolated to the faces.

    Fpre is checked when read: it must be dimensionless, calculated on every
    patch but the constraint ones, and have a positive determinant in every
    cell, and on every face once interpolated there. It is written at every
    write time at full precision, so that a continued run reads back the Fpre
    it was started with. checkLaws refuses the laws that cannot be
    prestressed, and is called by the manager before Fpre is read.

    Also declared here is finiteStrainPrestressedFields, the finite-strain
    kinematics a constitutive law sees when the solid is prestressed: the
    solver's deformation gradient composed with Fpre,

        F    = Fstep & Fpre
        Finv = inv(Fpre) & FinvStep
        J    = Jstep*det(Fpre)

    and the same for the old-time values, against the same Fpre. Fpre is on
    the right because it acts first: it takes the stress-free configuration
    to the mesh, and the solver's deformation gradient takes the mesh to the
    current configuration. The two do not commute once the principal
    directions differ between them.

    Only the laws see the composed values. The solver keeps its own,
    mesh-relative F and J for the momentum equation and its Piola transform,
    which is consistent because the laws return the Cauchy stress, and the
    Cauchy stress belongs to the current configuration whichever reference it
    was computed from.

    The composed fields are copies of the solver's, so they keep its patch
    types, and with them which patches are coupled; only the values change.
    A copy brings the old-time fields of what it copies with it, so the
    manager makes one of these once and recomposes it with update() at every
    later evaluation, rather than copying the solver's fields each time.

SourceFiles
    mechanicalConstitutiveLawPrestress.C

Author
    Conor Moore, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawPrestress_H
#define mechanicalConstitutiveLawPrestress_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "autoPtr.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawKinematicsFields.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class finiteStrainPrestressedFields Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class GeoMesh>
class finiteStrainPrestressedFields
{
public:

    typedef finiteStrainGeometricFields<PatchField, GeoMesh> geometricFields;
    typedef typename geometricFields::tensorFieldType tensorFieldType;
    typedef typename geometricFields::scalarFieldType scalarFieldType;


private:

    // Private Data

        // Declared before fields below, which refers to them

        //- Deformation gradient, Fstep & Fpre
        tensorFieldType F_;

        //- Old-time deformation gradient, F0step & Fpre
        tensorFieldType F0_;

        //- Inverse deformation gradient, inv(Fpre) & FinvStep
        tensorFieldType Finv_;

        //- Old-time inverse deformation gradient, inv(Fpre) & Finv0Step
        tensorFieldType Finv0_;

        //- Jacobian, Jstep*det(Fpre)
        scalarFieldType J_;

        //- Old-time Jacobian, J0step*det(Fpre)
        scalarFieldType J0_;


    // Private Member Functions

        //- An unregistered copy of fld, patch types and all
        template<class GeoField>
        static GeoField copyOf(const GeoField& fld)
        {
            return GeoField
            (
                IOobject
                (
                    fld.name() + "Prestressed",
                    fld.instance(),
                    fld.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fld
            );
        }

        //- Set out to a & b, internal and boundary values alike
        static void compose
        (
            tensorFieldType& out,
            const tensorFieldType& a,
            const tensorFieldType& b
        )
        {
            Foam::primitiveFieldRef(out) =
                Foam::primitiveField(a) & Foam::primitiveField(b);

            forAll(out.boundaryField(), patchI)
            {
                Foam::boundaryFieldRef(out)[patchI] ==
                (
                    a.boundaryField()[patchI] & b.boundaryField()[patchI]
                );
            }
        }

        //- Set out to a*b, internal and boundary values alike
        static void compose
        (
            scalarFieldType& out,
            const scalarFieldType& a,
            const scalarFieldType& b
        )
        {
            Foam::primitiveFieldRef(out) =
                Foam::primitiveField(a)*Foam::primitiveField(b);

            forAll(out.boundaryField(), patchI)
            {
                Foam::boundaryFieldRef(out)[patchI] ==
                (
                    a.boundaryField()[patchI]*b.boundaryField()[patchI]
                );
            }
        }

        //- No copy construct
        finiteStrainPrestressedFields
        (
            const finiteStrainPrestressedFields&
        ) = delete;

        //- No copy assignment
        void operator=(const finiteStrainPrestressedFields&) = delete;


public:

    //- The composed kinematics, on the fields above
    const geometricFields fields;


    // Constructors

        //- Compose the solver's kinematics, step, with the prestress
        //  deformation gradient Fpre, given with its inverse and determinant
        finiteStrainPrestressedFields
        (
            const geometricFields& step,
            const tensorFieldType& Fpre,
            const tensorFieldType& FpreInv,
            const scalarFieldType& Jpre
        )
        :
            F_(copyOf(step.F)),
            F0_(copyOf(step.F0)),
            Finv_(copyOf(step.Finv)),
            Finv0_(copyOf(step.Finv0)),
            J_(copyOf(step.J)),
            J0_(copyOf(step.J0)),
            fields(F_, F0_, Finv_, Finv0_, J_, J0_)
        {
            update(step, Fpre, FpreInv, Jpre);
        }


    // Member Functions

        //- Compose step with Fpre again, into the fields already made.
        //  step must be on the mesh, and have the patches, of the step this
        //  was constructed from; only its values are read
        void update
        (
            const geometricFields& step,
            const tensorFieldType& Fpre,
            const tensorFieldType& FpreInv,
            const scalarFieldType& Jpre
        )
        {
            compose(F_, step.F, Fpre);

            // The old-time values are composed with the same Fpre, which is
            // right for any law, since Fpre is constant during a run. A law
            // that builds its stress from increments, F & Finv0 say, would
            // cancel Fpre out; such laws are refused through admitsPrestress()
            compose(F0_, step.F0, Fpre);
            compose(Finv_, FpreInv, step.Finv);
            compose(Finv0_, FpreInv, step.Finv0);
            compose(J_, step.J, Jpre);
            compose(J0_, step.J0, Jpre);
        }
};


/*---------------------------------------------------------------------------*\
              Class mechanicalConstitutiveLawPrestress Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawPrestress
{
    // Private Classes

        //- Fpre, written at full round-trip precision (17 significant
        //  digits) whatever writePrecision the case asks for. This follows
        //  the rule set by the constitutive state, whose
        //  stateIOFieldProxy::writeData in
        //  mechanicalConstitutiveLawStateIOTemplates.C does the same:
        //  anything a restart reads back is written at full precision. A
        //  continued run reads Fpre back from the time it continues from,
        //  and a rounding of it would move the stress the laws see. writeData
        //  is overridden rather than writeObject, whose signature differs
        //  between the three forks
        class fullPrecisionVolTensorField
        :
            public volTensorField
        {
        public:

            fullPrecisionVolTensorField(const IOobject& io, const fvMesh& mesh)
            :
                volTensorField(io, mesh)
            {}

            virtual bool writeData(Ostream& os) const override
            {
                const int oldPrecision = os.precision(17);

                const bool ok = volTensorField::writeData(os);

                os.precision(oldPrecision);

                return ok;
            }
        };


    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Prestress deformation gradient, read from the start time, from
        //  the stress-free configuration to the mesh
        fullPrecisionVolTensorField Fpre_;

        //- inv(Fpre), taken once since Fpre does not change during a run
        volTensorField FpreInv_;

        //- det(Fpre), taken once since Fpre does not change during a run
        volScalarField Jpre_;

        //- Face Fpre, interpolated from the cells, made on demand
        mutable autoPtr<surfaceTensorField> FpreSurfacePtr_;

        //- Face inv(Fpre), the inverse of the face Fpre
        mutable autoPtr<surfaceTensorField> FpreInvSurfacePtr_;

        //- Face det(Fpre), the determinant of the face Fpre
        mutable autoPtr<surfaceScalarField> JpreSurfacePtr_;

        //- Time index at which the face fields were last made, so that on a
        //  moving mesh they are made again at most once per time step
        mutable label surfaceTimeIndex_;


    // Private Member Functions

        //- Check that Fpre is dimensionless and calculated on every patch
        //  but the constraint ones
        void check() const;

        //- Make the face fields, rebuilt when the mesh has moved, and check
        //  det(Fpre) on every face
        void makeSurfaceFields() const;

        //- The minimum of fld over its internal and boundary values on every
        //  processor
        template<class GeoField>
        static scalar gMinAll(const GeoField& fld);

        //- No copy construct
        mechanicalConstitutiveLawPrestress
        (
            const mechanicalConstitutiveLawPrestress&
        ) = delete;

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawPrestress&) = delete;


public:

    // Constructors

        //- Construct from the mesh, reading Fpre from the start time
        explicit mechanicalConstitutiveLawPrestress(const fvMesh& mesh);


    // Member Functions

        //- Refuse the laws that cannot be prestressed, and warn about the
        //  directions of those that take them in the stress-free
        //  configuration. Called before Fpre is read, so that a case that
        //  cannot be prestressed says so whatever is in the start time
        static void checkLaws
        (
            const PtrList<mechanicalConstitutiveLaw>& laws,
            const wordList& lawNames
        );

        //- Prestress deformation gradient
        const volTensorField& Fpre() const
        {
            return Fpre_;
        }

        //- inv(Fpre)
        const volTensorField& FpreInv() const
        {
            return FpreInv_;
        }

        //- det(Fpre)
        const volScalarField& Jpre() const
        {
            return Jpre_;
        }

        //- Face Fpre
        const surfaceTensorField& FpreSurface() const;

        //- Face inv(Fpre)
        const surfaceTensorField& FpreInvSurface() const;

        //- Face det(Fpre)
        const surfaceScalarField& JpreSurface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
