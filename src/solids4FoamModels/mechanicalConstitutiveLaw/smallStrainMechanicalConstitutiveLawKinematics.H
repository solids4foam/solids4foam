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

        //- Time step
        const scalar dt_;


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
            const UIndirectList<tensor>& gradD0,
            const scalar dt
        )
        :
            gradD_(gradD),
            gradD0_(gradD0),
            dt_(dt)
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

        //- Time increment
        scalar dt() const
        {
            return dt_;
        }

    // Edit

        //- Small strain tensor: epsilon = symm(gradD)
        inline symmTensor epsilon(const label i) const
        {
            return symm(gradD_[i]);
        }

        //- Strain increment: delta epsilon
        inline symmTensor deltaEpsilon(const label i) const
        {
            return symm(gradD_[i] - gradD0_[i]);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
