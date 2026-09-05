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
    Foam::GuccioneElasticMechanicalConstitutiveLaw

Description
    Guccione transversely isotropic hyperelastic law, for passive myocardium.

    The strain energy is exponential in a quadratic form of the Green-Lagrange
    strain, taken about a fibre direction:

        W = 0.5*k*(exp(Q) - 1)

        Q = ct*I1^2 - 2*ct*I2 + (cf - 2*cfs + ct)*I4^2 + 2*(cfs - ct)*I5

    where I1 and I2 are the first two invariants of the Green-Lagrange strain
    and I4 and I5 are the fibre invariants, formed with the structure tensor
    f0 (x) f0. Setting cf = cfs = ct makes Q isotropic and the fibre direction
    stops mattering.

    The fibre direction is prescribed state: a case supplies an f0 field, as
    setFibreField writes, and a case whose fibres are the same everywhere may
    give a single vector in this dictionary instead. It is normalised here, so
    it need not be given normalised, and a zero-length fibre is refused rather
    than divided by.

    Two things the legacy GuccioneElastic does are deliberately not here.

    The mixed displacement-pressure formulation is not ported. In the legacy
    law it is selected by a per-material dictionary entry, which makes a
    material's dictionary choose which equations the solver solves, and its two
    branches do not agree on the stress measure they return. That is being
    settled separately, and until it is this law is displacement only.

    Neither is the local coordinate system variant, which doubles the
    evaluation path for a case that no ported tutorial uses.

SourceFiles
    GuccioneElasticMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef GuccioneElasticMechanicalConstitutiveLaw_H
#define GuccioneElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class GuccioneElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class GuccioneElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- Density
        const dimensionedScalar rho_;

        //- Stiffness scale of the exponential
        const dimensionedScalar k_;

        //- Fibre, transverse and fibre-sheet coefficients of Q
        const scalar cf_;
        const scalar ct_;
        const scalar cfs_;

        //- Bulk modulus, a penalty enforcing near incompressibility
        const dimensionedScalar bulkModulus_;

        //- Small-strain shear modulus implied by k and the coefficients.
        //  Used for the implicit stiffness and nothing else
        const dimensionedScalar mu_;

        //- The fibre direction to use where a case supplies no f0 field.
        //  Zero when the case is expected to supply one, so that its absence
        //  is a zero-length fibre and refused rather than quietly isotropic
        const vector f0Default_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        GuccioneElasticMechanicalConstitutiveLaw
        (
            const GuccioneElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=
        (
            const GuccioneElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("GuccioneElastic");

    // Constructors

        //- Construct from dictionary
        explicit GuccioneElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~GuccioneElasticMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- Return the density
        virtual dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Return the bulk modulus
        virtual dimensionedScalar kappa() const override
        {
            return bulkModulus_;
        }

        //- Declare the fibre direction this law is given
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

        //- Finite-strain constitutive evaluation
        virtual void evaluate
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
