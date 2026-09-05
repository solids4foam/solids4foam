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

#include "anisotropicBiotElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "labelVector.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(anisotropicBiotElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        anisotropicBiotElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicBiotElasticMechanicalConstitutiveLaw::
anisotropicBiotElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    model2d_(false),
    A11_(0.0),
    A22_(0.0),
    A33_(0.0),
    A44_(0.0),
    A55_(0.0),
    A66_(0.0),
    A12_(0.0),
    A21_(0.0),
    A31_(0.0),
    A23_(0.0)
{
    // The directions the mesh solves in, injected by the manager because a
    // law is constructed from a dictionary and cannot ask the mesh
    const labelVector solutionD(dict.lookup("solutionD"));

    // A mesh is two-dimensional when z is the empty direction. The legacy law
    // tested this the other way round, so it ran the reduced plane model on
    // three-dimensional meshes and quietly ignored the out-of-plane constants
    // they supplied. Fixed in both places together, so the two still agree
    model2d_ = (solutionD[vector::Z] < 0);

    if (model2d_)
    {
        if (solutionD[vector::X] < 0 || solutionD[vector::Y] < 0)
        {
            FatalIOErrorInFunction(dict)
                << "For 2-D models, z must be the empty direction"
                << exit(FatalIOError);
        }

        // The reduced constants below are the plane stress ones: they come
        // from eliminating a zero out-of-plane stress, not a zero out-of-plane
        // strain. This law has no plane strain form, so a case that asked for
        // one is told rather than quietly given the other
        if (!dict.lookupOrDefault<Switch>("planeStress", false))
        {
            FatalIOErrorInFunction(dict)
                << "This law's two-dimensional reduction is plane stress, but "
                << "planeStress is not set in mechanicalProperties." << nl
                << "Set 'planeStress yes;', or use a law that offers a plane "
                << "strain reduction."
                << exit(FatalIOError);
        }


        const scalar Ex = readScalar(dict.lookup("Ex"));
        const scalar Ey = readScalar(dict.lookup("Ey"));
        const scalar vxy = readScalar(dict.lookup("nuxy"));
        const scalar Gxy = readScalar(dict.lookup("Gxy"));

        // Material constraints
        const scalar vyx = vxy*Ey/Ex;

        const scalar J = 1/(1 - vxy*vyx);
        A11_ = J*Ex;
        A22_ = J*Ey;
        A12_ = J*vyx*Ex;
        A21_ = J*vxy*Ey;
        A44_ = 2*Gxy;

        Info<< "    Orthotropic 2D stiffness: A11 " << A11_
            << ", A22 " << A22_ << ", A12 " << A12_
            << ", A21 " << A21_ << ", A44 " << A44_ << endl;
    }
    else
    {
        const scalar Ex = readScalar(dict.lookup("Ex"));
        const scalar Ey = readScalar(dict.lookup("Ey"));
        const scalar Ez = readScalar(dict.lookup("Ez"));
        const scalar vxy = readScalar(dict.lookup("nuxy"));
        const scalar vyz = readScalar(dict.lookup("nuyz"));
        const scalar vzx = readScalar(dict.lookup("nuzx"));
        const scalar Gxy = readScalar(dict.lookup("Gxy"));
        const scalar Gyz = readScalar(dict.lookup("Gyz"));
        const scalar Gzx = readScalar(dict.lookup("Gzx"));

        // Material constraints
        const scalar vyx = vxy*Ey/Ex;
        const scalar vxz = vzx*Ex/Ez;
        const scalar vzy = vyz*Ez/Ey;

        const scalar J =
            (1.0 - vxy*vyx - vyz*vzy - vzx*vxz - 2*vyx*vzy*vxz)/(Ex*Ey*Ez);
        A11_ = (1.0 - vyz*vzy)/(J*Ey*Ez);
        A22_ = (1.0 - vxz*vzx)/(J*Ex*Ez);
        A33_ = (1.0 - vyx*vxy)/(J*Ey*Ex);
        A12_ = (vxy + vzy*vxz)/(J*Ex*Ez);
        A31_ = (vzx + vyx*vzy)/(J*Ey*Ez);
        A23_ = (vyz + vyx*vxz)/(J*Ex*Ey);
        A44_ = 2*Gxy;
        A55_ = 2*Gyz;
        A66_ = 2*Gzx;

        Info<< "    Orthotropic 3D stiffness: A11 " << A11_
            << ", A22 " << A22_ << ", A33 " << A33_
            << ", A12 " << A12_ << ", A31 " << A31_
            << ", A23 " << A23_ << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::anisotropicBiotElasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const UIndirectList<tensor>& gradD = kin.gradD();

    forAll(sigma, i)
    {
        const symmTensor epsilon(symm(gradD[i]));

        const scalar e11 = epsilon[symmTensor::XX];
        const scalar e22 = epsilon[symmTensor::YY];
        const scalar e33 = epsilon[symmTensor::ZZ];
        const scalar e12 = epsilon[symmTensor::XY];
        const scalar e23 = epsilon[symmTensor::YZ];
        const scalar e31 = epsilon[symmTensor::XZ];

        if (model2d_)
        {
            sigma[i][symmTensor::XX] = A11_*e11 + A12_*e22;
            sigma[i][symmTensor::YY] = A21_*e11 + A22_*e22;
            sigma[i][symmTensor::XY] = A44_*e12;

            // The reduced constants above are the plane stress ones, so the
            // out-of-plane components are zero by construction. The legacy law
            // left them untouched instead, which came to the same thing only
            // because its stress field starts at zero and nothing else writes
            // them; here it is stated, so the result does not depend on what
            // storage the law was handed
            sigma[i][symmTensor::ZZ] = 0.0;
            sigma[i][symmTensor::YZ] = 0.0;
            sigma[i][symmTensor::XZ] = 0.0;
        }
        else
        {
            sigma[i][symmTensor::XX] = A11_*e11 + A12_*e22 + A31_*e33;
            sigma[i][symmTensor::YY] = A12_*e11 + A22_*e22 + A23_*e33;
            sigma[i][symmTensor::ZZ] = A31_*e11 + A23_*e22 + A33_*e33;
            sigma[i][symmTensor::XY] = A44_*e12;
            sigma[i][symmTensor::YZ] = A55_*e23;
            sigma[i][symmTensor::XZ] = A66_*e31;
        }
    }

    // Scalar tangent: only if explicitly requested. The legacy impK uses the
    // largest direct stiffness
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        const scalar Keff = max(A11_, max(A22_, A33_));

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }
}


// ************************************************************************* //
