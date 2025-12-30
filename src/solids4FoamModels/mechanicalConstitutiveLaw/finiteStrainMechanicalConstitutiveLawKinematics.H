/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::finiteStrainMechanicalConstitutiveLawKinematics

Description
    Lightweight container class holding the kinematic fields required by
    finite-strain mechanical constitutive laws.

    The class provides non-owning views (via UIndirectList) onto stored
    finite-strain kinematic quantities, including the deformation gradient,
    its inverse, and its determinant, at the current and previous time
    levels, together with the time increment.

    All kinematic quantities are assumed to be computed, updated, and
    stored externally by the solidModel and mechanicalConstitutiveLawManager.
    This class does not perform any kinematic updates or time integration.

    No strain or stress measures are stored explicitly. Instead, any derived
    quantities required by a constitutive law (e.g. relative deformation
    gradients or strain measures) are computed on demand from the stored
    kinematic fields.

    This class is intentionally passive and allocation-free. Its sole
    responsibility is to present the finite-strain kinematic state required
    by mechanical constitutive laws in a clear, consistent, and efficient
    manner.

    Usage:
    An instance of finiteStrainMechanicalConstitutiveLawKinematics is
    constructed by the mechanicalConstitutiveLawManager using UIndirectList
    views of the finite-strain kinematic fields and passed to a
    mechanicalConstitutiveLaw during constitutive evaluation.

    Notes:
    - The class is non-owning and relies on the lifetime of the referenced
      fields exceeding that of the kinematics object.
    - The separation between stored kinematic measures and derived
      quantities avoids duplication of data and allows constitutive laws
      to implement a wide range of formulations (e.g. total or updated
      Lagrangian, hyperelasticity, finite-strain plasticity).
    - This class is specific to finite-strain formulations; small-strain
      formulations are handled by a separate kinematics container.

SourceFiles
    finiteStrainMechanicalConstitutiveLawKinematics.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef finiteStrainMechanicalConstitutiveLawKinematics_H
#define finiteStrainMechanicalConstitutiveLawKinematics_H

#include "finiteStrainMechanicalConstitutiveLawKinematics.H"
#include "UIndirectList.H"
#include "tensor.H"
#include "symmTensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class finiteStrainMechanicalConstitutiveLawKinematics Declaration
\*---------------------------------------------------------------------------*/

class finiteStrainMechanicalConstitutiveLawKinematics
{
    // Private Data

        //- Deformation gradient at current time
        const UIndirectList<tensor>& F_;

        //- Deformation gradient at previous time
        const UIndirectList<tensor>& F0_;

        //- Determinant of deformation gradient
        const UIndirectList<scalar>& J_;

        //- Determinant of deformation gradient at previous time
        const UIndirectList<scalar>& J0_;

        //- Inverse deformation gradient
        const UIndirectList<tensor>& Finv_;

        //- Inverse deformation gradient at previous time
        const UIndirectList<tensor>& Finv0_;

        //- Time step
        const scalar dt_;


    // Private Member Functions

        //- No copy construct
        finiteStrainMechanicalConstitutiveLawKinematics
        (
            const finiteStrainMechanicalConstitutiveLawKinematics&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const finiteStrainMechanicalConstitutiveLawKinematics&
        ) = delete;


public:


    // Constructors

        //- Construct from current and old displacement gradients
        finiteStrainMechanicalConstitutiveLawKinematics
        (
            const UIndirectList<tensor>& F,
            const UIndirectList<tensor>& F0,
            const UIndirectList<scalar>& J,
            const UIndirectList<scalar>& J0,
            const UIndirectList<tensor>& Finv,
            const UIndirectList<tensor>& Finv0,
            const scalar dt
        )
        :
            F_(F),
            F0_(F0),
            J_(J),
            J0_(J0),
            Finv_(Finv),
            Finv0_(Finv0),
            dt_(dt)
        {}


    //- Destructor
    ~finiteStrainMechanicalConstitutiveLawKinematics()
    {}


    // Member Functions

    // Access

        //- Access current displacement gradient
        const UIndirectList<tensor>& F() const
        {
            return F_;
        }

        //- Access previous time displacement gradient
        const UIndirectList<tensor>& F0() const
        {
            return F0_;
        }

        //- Access current displacement gradient determinant
        const UIndirectList<scalar>& J() const
        {
            return J_;
        }

        //- Access previous time displacement gradient determinant
        const UIndirectList<scalar>& J0() const
        {
            return J0_;
        }

        //- Access current displacement gradient inverse
        const UIndirectList<tensor>& Finv() const
        {
            return Finv_;
        }

        //- Access previous time displacement gradient inverse
        const UIndirectList<tensor>& Finv0() const
        {
            return Finv0_;
        }

        //- Time increment
        scalar dt() const
        {
            return dt_;
        }

    // Edit

        //- Relative deformation gradient: F*inv(F0)
        inline tensor relativeF(const label i) const
        {
            return F_[i] & inv(F0_[i]);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
