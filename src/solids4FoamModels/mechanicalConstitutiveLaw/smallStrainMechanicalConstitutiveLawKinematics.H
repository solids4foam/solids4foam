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
    Foam::smallStrainMechanicalConstitutiveLawKinematics

Description
    Lightweight container class holding the kinematic fields required by
    small-strain mechanical constitutive laws.

    The class provides non-owning views (via UIndirectList) onto the current
    and previous displacement gradient fields, together with the time
    increment. These quantities are assumed to be computed and managed
    externally by the solidModel and mechanicalConstitutiveLawManager.

    No strain fields are stored explicitly. Instead, small-strain
    measures such as the symmetric strain tensor and strain increment
    are provided as derived quantities computed on demand from the
    displacement gradient.

    This class is intentionally passive and performs no kinematic
    updates or time integration. Its sole responsibility is to present
    the kinematic state required by small-strain mechanical constitutive
    laws in a clear, consistent, and allocation-free manner.

    Usage:
    An instance of smallStrainMechanicalConstitutiveLawKinematics is
    constructed by the mechanicalConstitutiveLawManager using UIndirectList
    views of the displacement gradient fields and passed to a
    mechanicalConstitutiveLaw during constitutive evaluation.

    Notes:
    - The class is non-owning and relies on the lifetime of the referenced
      fields exceeding that of the kinematics object.
    - The separation between stored displacement gradients and derived
      strain measures avoids duplication of kinematic storage and allows
      constitutive laws to choose the most appropriate formulation
      (total, incremental, or rate-based).
    - This class is specific to small-strain formulations; finite-strain
      formulations are handled by a separate kinematics container.

SourceFiles
    smallStrainMechanicalConstitutiveLawKinematics.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef smallStrainMechanicalConstitutiveLawKinematics_H
#define smallStrainMechanicalConstitutiveLawKinematics_H

#include "smallStrainMechanicalConstitutiveLawKinematics.H"
#include "UIndirectList.H"
#include "tensor.H"
#include "symmTensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class smallStrainMechanicalConstitutiveLawKinematics Declaration
\*---------------------------------------------------------------------------*/

class smallStrainMechanicalConstitutiveLawKinematics
{
    // Private Data

        //- Displacement gradient at current time
        const UIndirectList<tensor>& gradD_;

        //- Displacement gradient at previous time
        //  This may be unused for some laws, e.g. purely elastic laws
        const UIndirectList<tensor>& gradD0_;



    // Private Member Functions

        //- No copy construct
        smallStrainMechanicalConstitutiveLawKinematics
        (
            const smallStrainMechanicalConstitutiveLawKinematics&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const smallStrainMechanicalConstitutiveLawKinematics&
        ) = delete;


public:


    // Constructors

        //- Construct from current and old displacement gradients
        smallStrainMechanicalConstitutiveLawKinematics
        (
            const UIndirectList<tensor>& gradD,
            const UIndirectList<tensor>& gradD0
        )
        :
            gradD_(gradD),
            gradD0_(gradD0)
        {}


    //- Destructor
    ~smallStrainMechanicalConstitutiveLawKinematics()
    {}


    // Member Functions

    // Access

        //- Access current displacement gradient
        const UIndirectList<tensor>& gradD() const
        {
            return gradD_;
        }

        //- Access old displacement gradient
        const UIndirectList<tensor>& gradD0() const
        {
            return gradD0_;
        }


    // Edit

        //- Return a perturbed copy of gradD at ipI
        //  Voigt strain perturbation to integration point ipI
        //  j = Voigt index [0..5]
        //  h = perturbation magnitude
        //  This can be used for calculating material tangents with finite
        //  differences
        //- Return gradD at integration point ipI with Voigt strain component
        //  j perturbed by h.
        //  The perturbation is defined on the Voigt strain vector, which uses
        //  engineering shear: gamma_xy = 2*eps_xy. A shear component therefore
        //  moves each off-diagonal of gradD by h/2, so that the corresponding
        //  Voigt component moves by exactly h. This is what makes a
        //  finite-difference tangent directly comparable with an analytical
        //  one, where for isotropic elasticity C(XY,XY) is mu rather than 2*mu
        inline tensor gradDPerturbed
        (
            const label ipI, const label j, const scalar h
        ) const
        {
            #ifdef FULLDEBUG
            if (j < 0 || j >= 6)
            {
                FatalErrorInFunction
                    << "Invalid Voigt index j = " << j
                    << " (expected 0..5)"
                    << exit(FatalError);
            }
            #endif

            // Copy of gradD
            tensor g = gradD_[ipI];

            switch (j)
            {
                case symmTensor::XX:
                    g.xx() += h;
                    break;

                case symmTensor::YY:
                    g.yy() += h;
                    break;

                case symmTensor::ZZ:
                    g.zz() += h;
                    break;

                case symmTensor::YZ:
                    g.yz() += 0.5*h;
                    g.zy() += 0.5*h;
                    break;

                case symmTensor::XZ:
                    g.xz() += 0.5*h;
                    g.zx() += 0.5*h;
                    break;

                case symmTensor::XY:
                    g.xy() += 0.5*h;
                    g.yx() += 0.5*h;
                    break;

                default:
                    FatalErrorInFunction
                        << "Invalid Voigt component index j = " << j
                            << exit(FatalError);
            }

            return g;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
