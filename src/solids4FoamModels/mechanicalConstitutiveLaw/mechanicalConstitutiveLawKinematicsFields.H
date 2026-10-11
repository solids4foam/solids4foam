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
    Foam::smallStrainKinematicsFields, Foam::smallStrainKinematicsViews,
    Foam::finiteStrainKinematicsFields, Foam::finiteStrainKinematicsViews,
    Foam::smallStrainGeometricFields, Foam::finiteStrainGeometricFields

Description
    The kinematic inputs of an evaluation, as lists over integration points,
    and the views onto them that a law's kinematics wrapper is built from.

    The mechanicalConstitutiveLawManager evaluates small- and finite-strain
    laws the same way; what differs is which lists the kinematics carry and
    which convergence scale a law reports for them. A Fields struct holds the
    lists and says which scale; its Views class addresses them at one law's
    points and holds the kinematics built on those views. The kinematics keep
    references to the views and cannot be copied, so a Views object is built
    in place and never returned by value; its members are declared in the
    order they must be constructed in.

    kinematicsViewsOf<Fields>::type names the Views class of a Fields, for
    code written once for both.

    The GeometricFields structs hold the same inputs as volFields or
    surfaceFields, so the manager's field overloads can be written once for
    both strain measures. Each gives its internal and patch values as a
    Fields, checks its fields' meshes, and names the overloads it serves.

SourceFiles
    (header only)

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawKinematicsFields_H
#define mechanicalConstitutiveLawKinematicsFields_H

#include "smallStrainMechanicalConstitutiveLawKinematics.H"
#include "finiteStrainMechanicalConstitutiveLawKinematics.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawState.H"
#include "UIndirectList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Struct smallStrainKinematicsFields Declaration
\*---------------------------------------------------------------------------*/

struct smallStrainKinematicsFields
{
    // The lists are held by reference, as are the views built on them, so
    // each must outlive both: pass a field that exists, never a temporary
    // such as a tmp or a patchInternalField()

    typedef smallStrainMechanicalConstitutiveLawKinematics kinematicsType;

    const UList<tensor>& gradD;
    const UList<tensor>& gradD0;

    smallStrainKinematicsFields
    (
        const UList<tensor>& gradDIn,
        const UList<tensor>& gradD0In
    )
    :
        gradD(gradDIn),
        gradD0(gradD0In)
    {}

    //- The scale a law normalises its convergence test by, over these points
    static scalar convergenceScale
    (
        const mechanicalConstitutiveLaw& law,
        const kinematicsType& kin,
        const mechanicalConstitutiveLawState& state
    )
    {
        return law.smallStrainConvergenceScale(kin, state);
    }

    //- The name the flat-list evaluation reports itself by
    static const char* flatContext()
    {
        return "updateStressSmallStrain (flat list)";
    }

    //- What the flat-list evaluation suggests when it cannot collapse a
    //  shared point: the face and point overloads, which can
    static const char* multiLawHint()
    {
        return "\nUse the surfaceField or pointField overload, which takes a "
            "stressCollapseRule.";
    }
};


/*---------------------------------------------------------------------------*\
                   Class smallStrainKinematicsViews Declaration
\*---------------------------------------------------------------------------*/

class smallStrainKinematicsViews
{
    const UIndirectList<tensor> gradD_;
    const UIndirectList<tensor> gradD0_;

public:

    //- The kinematics on the views above, which it refers to
    const smallStrainMechanicalConstitutiveLawKinematics kin;

    smallStrainKinematicsViews
    (
        const smallStrainKinematicsFields& fields,
        const labelUList& ids
    )
    :
        gradD_(fields.gradD, ids),
        gradD0_(fields.gradD0, ids),
        kin(gradD_, gradD0_)
    {}
};


/*---------------------------------------------------------------------------*\
                  Struct finiteStrainKinematicsFields Declaration
\*---------------------------------------------------------------------------*/

struct finiteStrainKinematicsFields
{
    // The lists are held by reference, as are the views built on them, so
    // each must outlive both: pass a field that exists, never a temporary
    // such as a tmp or a patchInternalField()

    typedef finiteStrainMechanicalConstitutiveLawKinematics kinematicsType;

    const UList<tensor>& F;
    const UList<tensor>& F0;
    const UList<tensor>& Finv;
    const UList<tensor>& Finv0;
    const UList<scalar>& J;
    const UList<scalar>& J0;

    finiteStrainKinematicsFields
    (
        const UList<tensor>& FIn,
        const UList<tensor>& F0In,
        const UList<tensor>& FinvIn,
        const UList<tensor>& Finv0In,
        const UList<scalar>& JIn,
        const UList<scalar>& J0In
    )
    :
        F(FIn),
        F0(F0In),
        Finv(FinvIn),
        Finv0(Finv0In),
        J(JIn),
        J0(J0In)
    {}

    //- The scale a law normalises its convergence test by, over these points
    static scalar convergenceScale
    (
        const mechanicalConstitutiveLaw& law,
        const kinematicsType& kin,
        const mechanicalConstitutiveLawState& state
    )
    {
        return law.localConvergenceScale(kin, state);
    }

    //- The name the flat-list evaluation reports itself by
    static const char* flatContext()
    {
        return "updateStressFiniteStrain (flat list)";
    }

    //- There is no finite-strain overload that collapses a shared point, so
    //  nothing to suggest
    static const char* multiLawHint()
    {
        return "";
    }
};


/*---------------------------------------------------------------------------*\
                  Class finiteStrainKinematicsViews Declaration
\*---------------------------------------------------------------------------*/

class finiteStrainKinematicsViews
{
    const UIndirectList<tensor> F_;
    const UIndirectList<tensor> F0_;
    const UIndirectList<tensor> Finv_;
    const UIndirectList<tensor> Finv0_;
    const UIndirectList<scalar> J_;
    const UIndirectList<scalar> J0_;

public:

    //- The kinematics on the views above, which it refers to
    const finiteStrainMechanicalConstitutiveLawKinematics kin;

    finiteStrainKinematicsViews
    (
        const finiteStrainKinematicsFields& fields,
        const labelUList& ids
    )
    :
        F_(fields.F, ids),
        F0_(fields.F0, ids),
        Finv_(fields.Finv, ids),
        Finv0_(fields.Finv0, ids),
        J_(fields.J, ids),
        J0_(fields.J0, ids),
        kin(F_, F0_, Finv_, Finv0_, J_, J0_)
    {}
};


/*---------------------------------------------------------------------------*\
                      Struct kinematicsViewsOf Declaration
\*---------------------------------------------------------------------------*/

template<class Fields>
struct kinematicsViewsOf;

template<>
struct kinematicsViewsOf<smallStrainKinematicsFields>
{
    typedef smallStrainKinematicsViews type;
};

template<>
struct kinematicsViewsOf<finiteStrainKinematicsFields>
{
    typedef finiteStrainKinematicsViews type;
};


//- Check a field is defined on the expected mesh
template<class GeoField>
inline void checkFieldMesh(const polyMesh& mesh, const GeoField& field)
{
    if (&field.mesh() != &mesh)
    {
        FatalError
            << "Field '" << field.name() << "' is associated with a "
            << "different mesh." << nl
            << "Expected mesh: " << mesh.name() << nl
            << "Field mesh:    " << field.mesh().name()
            << exit(FatalError);
    }
}


/*---------------------------------------------------------------------------*\
                   Struct smallStrainGeometricFields Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class GeoMesh>
struct smallStrainGeometricFields
{
    // The fields are held by reference, as are the lists taken from them, so
    // each must outlive both

    typedef smallStrainKinematicsFields listFields;
    typedef GeometricField<tensor, PatchField, GeoMesh> tensorFieldType;

    const tensorFieldType& gradD;
    const tensorFieldType& gradD0;

    smallStrainGeometricFields
    (
        const tensorFieldType& gradDIn,
        const tensorFieldType& gradD0In
    )
    :
        gradD(gradDIn),
        gradD0(gradD0In)
    {}

    //- The name the stress update reports itself by
    static const char* stressContext()
    {
        return "updateStressSmallStrain";
    }

    //- The name the scalar tangent update reports itself by
    static const char* scalarTangentContext()
    {
        return "updateScalarTangent";
    }

    //- Check each field is defined on mesh
    void checkMeshes(const polyMesh& mesh) const
    {
        checkFieldMesh(mesh, gradD);
        checkFieldMesh(mesh, gradD0);
    }

    //- The field whose patches say which are coupled
    const tensorFieldType& lead() const
    {
        return gradD;
    }

    //- The internal values
    listFields internal() const
    {
        return listFields
        (
            Foam::primitiveField(gradD),
            Foam::primitiveField(gradD0)
        );
    }

    //- The values on patch patchI
    listFields patch(const label patchI) const
    {
        return listFields
        (
            gradD.boundaryField()[patchI],
            gradD0.boundaryField()[patchI]
        );
    }
};


/*---------------------------------------------------------------------------*\
                  Struct finiteStrainGeometricFields Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class GeoMesh>
struct finiteStrainGeometricFields
{
    // The fields are held by reference, as are the lists taken from them, so
    // each must outlive both

    typedef finiteStrainKinematicsFields listFields;
    typedef GeometricField<tensor, PatchField, GeoMesh> tensorFieldType;
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarFieldType;

    const tensorFieldType& F;
    const tensorFieldType& F0;
    const tensorFieldType& Finv;
    const tensorFieldType& Finv0;
    const scalarFieldType& J;
    const scalarFieldType& J0;

    finiteStrainGeometricFields
    (
        const tensorFieldType& FIn,
        const tensorFieldType& F0In,
        const tensorFieldType& FinvIn,
        const tensorFieldType& Finv0In,
        const scalarFieldType& JIn,
        const scalarFieldType& J0In
    )
    :
        F(FIn),
        F0(F0In),
        Finv(FinvIn),
        Finv0(Finv0In),
        J(JIn),
        J0(J0In)
    {}

    //- The name the stress update reports itself by
    static const char* stressContext()
    {
        return "updateStressFiniteStrain";
    }

    //- The name the scalar tangent update reports itself by
    static const char* scalarTangentContext()
    {
        return "updateScalarTangentFiniteStrain";
    }

    //- Check each field is defined on mesh
    void checkMeshes(const polyMesh& mesh) const
    {
        checkFieldMesh(mesh, F);
        checkFieldMesh(mesh, F0);
        checkFieldMesh(mesh, Finv);
        checkFieldMesh(mesh, Finv0);
        checkFieldMesh(mesh, J);
        checkFieldMesh(mesh, J0);
    }

    //- The field whose patches say which are coupled
    const tensorFieldType& lead() const
    {
        return F;
    }

    //- The internal values
    listFields internal() const
    {
        return listFields
        (
            Foam::primitiveField(F),
            Foam::primitiveField(F0),
            Foam::primitiveField(Finv),
            Foam::primitiveField(Finv0),
            Foam::primitiveField(J),
            Foam::primitiveField(J0)
        );
    }

    //- The values on patch patchI
    listFields patch(const label patchI) const
    {
        return listFields
        (
            F.boundaryField()[patchI],
            F0.boundaryField()[patchI],
            Finv.boundaryField()[patchI],
            Finv0.boundaryField()[patchI],
            J.boundaryField()[patchI],
            J0.boundaryField()[patchI]
        );
    }
};


typedef smallStrainGeometricFields<fvPatchField, volMesh>
    smallStrainVolFields;
typedef smallStrainGeometricFields<fvsPatchField, surfaceMesh>
    smallStrainSurfaceFields;
typedef finiteStrainGeometricFields<fvPatchField, volMesh>
    finiteStrainVolFields;
typedef finiteStrainGeometricFields<fvsPatchField, surfaceMesh>
    finiteStrainSurfaceFields;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
