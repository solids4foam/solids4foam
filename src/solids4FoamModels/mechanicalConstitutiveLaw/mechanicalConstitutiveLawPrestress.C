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

#include "mechanicalConstitutiveLawPrestress.H"
#include "calculatedFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::scalar Foam::mechanicalConstitutiveLawPrestress::gMinAll
(
    const GeoField& fld
)
{
    scalar minValue = gMin(Foam::primitiveField(fld));

    forAll(fld.boundaryField(), patchI)
    {
        const scalarField& pfld = fld.boundaryField()[patchI];

        if (pfld.size())
        {
            minValue = min(minValue, min(pfld));
        }
    }

    reduce(minValue, minOp<scalar>());

    return minValue;
}


void Foam::mechanicalConstitutiveLawPrestress::check() const
{
    if (Fpre_.dimensions() != dimless)
    {
        FatalErrorInFunction
            << "Fpre is a deformation gradient and must be dimensionless, "
            << "but it has dimensions " << Fpre_.dimensions()
            << exit(FatalError);
    }

    // The patch values of Fpre are composed with the solver's, as
    // updatePrestress writes them, and a patch type that evaluates its own
    // value, such as zeroGradient or fixedValue, would replace or freeze them.
    // So Fpre is calculated on every patch but the constraint ones (coupled,
    // empty, symmetry, wedge), whose values follow from the cells
    OStringStream refused;

    forAll(Fpre_.boundaryField(), patchI)
    {
        const fvPatchTensorField& pFpre = Fpre_.boundaryField()[patchI];

        if
        (
            !polyPatch::constraintType(pFpre.patch().type())
         && pFpre.type() != calculatedFvPatchField<tensor>::typeName
        )
        {
            refused
                << "    patch " << pFpre.patch().name() << " is "
                << pFpre.type() << nl;
        }
    }

    if (!refused.str().empty())
    {
        FatalErrorInFunction
            << "Fpre must be calculated on every patch but the "
            << "constraint ones, and here:" << nl
            << refused.str().c_str()
            << "A patch type that sets its own value would replace the "
            << "patch values of Fpre that updatePrestress composes"
            << exit(FatalError);
    }
}


void Foam::mechanicalConstitutiveLawPrestress::makeSurfaceFields() const
{
    // Made once on a static mesh. On a moving one the interpolation weights
    // change, so it is made again, and the faces checked again, once in each
    // time step
    if
    (
        FpreSurfacePtr_.valid()
     && (
            !mesh_.moving()
         || surfaceTimeIndex_ == mesh_.time().timeIndex()
        )
    )
    {
        return;
    }

    surfaceTimeIndex_ = mesh_.time().timeIndex();

    // The inverse and determinant are taken of the interpolated Fpre,
    // not interpolated themselves, so that the three agree with one another
    // on every face
    FpreSurfacePtr_.reset
    (
        new surfaceTensorField
        (
            IOobject
            (
                "Fpref",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvc::interpolate(Fpre_)
        )
    );

    JpreSurfacePtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Jpref",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            det(autoPtrRef(FpreSurfacePtr_))
        )
    );

    // A positive det(Fpre) in every cell does not make it positive on every
    // face: I and diag(-3, -1/3, 1) both have det 1, and their mean has
    // det -1/3. So the faces are checked, processor patches included,
    // before the face Fpre is inverted
    const scalar minJpref = gMinAll(autoPtrRef(JpreSurfacePtr_));

    if (minJpref <= SMALL)
    {
        FatalErrorInFunction
            << "det(Fpre) must be positive on every face, but the minimum "
            << "of det(Fpre) interpolated to the faces is " << minJpref
            << "." << nl
            << "Fpre is positive in every cell but changes too sharply "
            << "between neighbouring cells to be interpolated"
            << exit(FatalError);
    }

    FpreInvSurfacePtr_.reset
    (
        new surfaceTensorField
        (
            IOobject
            (
                "FpreInvf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            inv(autoPtrRef(FpreSurfacePtr_))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawPrestress::mechanicalConstitutiveLawPrestress
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    Fpre_
    (
        IOobject
        (
            "Fpre",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    FpreInv_
    (
        IOobject
        (
            "FpreInv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    Jpre_
    (
        IOobject
        (
            "Jpre",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    FpreSurfacePtr_(),
    FpreInvSurfacePtr_(),
    JpreSurfacePtr_(),
    surfaceTimeIndex_(-1)
{
    Info<< "Prestress: the mechanical constitutive laws see the deformation "
        << "gradient composed with the prestress deformation gradient Fpre, "
        << "read from time "
        << mesh_.time().timeName() << endl;

    check();

    // Fills the processor patches, which the file does not decide
    Fpre_.correctBoundaryConditions();

    Jpre_ = det(Fpre_);

    // Checked before the inverse is taken: an Fpre that turns an element
    // inside out has no inverse, and one that is singular has no stress-free
    // configuration to refer to
    const scalar minJpre = gMinAll(Jpre_);

    if (minJpre <= SMALL)
    {
        FatalErrorInFunction
            << "det(Fpre) must be positive everywhere, but its minimum is "
            << minJpre << exit(FatalError);
    }

    FpreInv_ = inv(Fpre_);

    Info<< "    det(Fpre) ranges from " << minJpre << " to "
        << gMax(Foam::primitiveField(Jpre_)) << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mechanicalConstitutiveLawPrestress::checkLaws
(
    const PtrList<mechanicalConstitutiveLaw>& laws,
    const wordList& lawNames
)
{
    OStringStream refused;

    forAll(laws, lawI)
    {
        if (!laws[lawI].admitsPrestress())
        {
            refused
                << "    material " << lawNames[lawI] << " uses "
                << laws[lawI].type() << nl;
        }
    }

    if (!refused.str().empty())
    {
        FatalErrorInFunction
            << "prestress is on, but these laws cannot be prestressed:"
            << nl << refused.str().c_str()
            << "Prestress applies to finite-strain laws whose stress is a "
            << "function of the total deformation gradient. A law that "
            << "integrates an increment from a history variable, such as "
            << "a plastic one, cancels Fpre out of the increment, and a "
            << "small-strain law has no deformation gradient to compose "
            << "it with" << exit(FatalError);
    }

    // A law with directions, such as fibres, takes them in the
    // stress-free configuration, and they are not mapped back from the
    // imaged geometry they were usually built on
    forAll(laws, lawI)
    {
        const wordList directions(laws[lawI].stressFreeDirections());

        if (directions.size())
        {
            // Written as "Ec Ea" rather than as the list "2(Ec Ea)"
            string names;

            forAll(directions, dirI)
            {
                names += (dirI ? " " : "") + directions[dirI];
            }

            WarningInFunction
                << "Material " << lawNames[lawI] << " uses "
                << laws[lawI].type() << ", which takes "
                << names.c_str() << " as directions in the stress-free "
                << "configuration." << nl
                << "    With prestress they are not mapped back from the "
                << "imaged geometry, so in the imaged geometry the law's "
                << "directions lie along Fpre & d rather than along d. "
                << "updatePrestress reports the angle between the two on "
                << "each pass" << nl << endl;
        }
    }
}


const Foam::surfaceTensorField&
Foam::mechanicalConstitutiveLawPrestress::FpreSurface() const
{
    makeSurfaceFields();

    return autoPtrRef(FpreSurfacePtr_);
}


const Foam::surfaceTensorField&
Foam::mechanicalConstitutiveLawPrestress::FpreInvSurface() const
{
    makeSurfaceFields();

    return autoPtrRef(FpreInvSurfacePtr_);
}


const Foam::surfaceScalarField&
Foam::mechanicalConstitutiveLawPrestress::JpreSurface() const
{
    makeSurfaceFields();

    return autoPtrRef(JpreSurfacePtr_);
}


// ************************************************************************* //
