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
    Foam::finiteStrainPrestressedFields

Description
    The finite-strain kinematics a constitutive law sees when the solid is
    prestressed: the solver's deformation gradient composed with a fixed
    prestress deformation gradient Fpre, which maps a stress-free
    configuration onto the mesh the solver is solving on.

        F    = Fstep & Fpre
        Finv = inv(Fpre) & Finv_step
        J    = J_step*det(Fpre)

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

SourceFiles
    (header only)

Author
    Conor Moore, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawPrestress_H
#define mechanicalConstitutiveLawPrestress_H

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

        tensorFieldType F_;
        tensorFieldType F0_;
        tensorFieldType Finv_;
        tensorFieldType Finv0_;
        scalarFieldType J_;
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
            compose(F_, step.F, Fpre);
            compose(F0_, step.F0, Fpre);
            compose(Finv_, FpreInv, step.Finv);
            compose(Finv0_, FpreInv, step.Finv0);
            compose(J_, step.J, Jpre);
            compose(J0_, step.J0, Jpre);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
