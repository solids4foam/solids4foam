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
    Foam::MooneyRivlinElasticMechanicalConstitutiveLaw

Description
    Finite-strain isotropic Mooney-Rivlin hyperelastic constitutive law.

    The isochoric response is

        isoB = J^(-2/3) F F^T,   I1 = tr(isoB),
        I2   = (I1^2 - tr(isoB.isoB))/2,

        s = 2 (c10 + c11 (I2 - 3)) isoB - 2 (c01 + c11 (I1 - 3)) isoB^-1,

    and the Cauchy stress is

        sigma = ( dev(s) + p I ) / J,   p = K (J^2 - 1)/2.

    The material is given as c10, c01 and c11, together with either the bulk
    modulus K or Poisson's ratio nu. The shear modulus follows from the pure
    shear state as mu = 2 (c10 + c01), and where nu is given the bulk modulus
    is taken from linear elasticity with E = 6 (c10 + c01).

    Note on the hydrostatic stress. The legacy MooneyRivlinElastic law passes
    its hydrostatic stress through updateSigmaHyd, which optionally solves a
    Laplacian smoothing equation over the whole field. That is deliberately
    NOT done here: a mechanicalConstitutiveLaw is a pure function of
    integration-point kinematics and old-time state, and a field-level
    implicit solve does not belong in that interface. Pressure smoothing
    stabilises the discretisation rather than describing the material, so it
    belongs to the solid model, where it is to be provided separately. This
    law therefore returns the unsmoothed p = K (J^2 - 1)/2.

    The legacy law also adds a pushed-forward initial stress,
    symm(F sigma0 F^T), to the Cauchy stress. No law in this framework carries
    an initial stress yet, so that term is absent here; a case that relies on
    a non-zero sigma0 is not yet reproduced.

    The material constants are uniform here, whereas the legacy law holds them
    as volScalarFields that may be read per cell. Spatially varying material
    constants are part of the wider question of how the framework carries
    per-integration-point material data, along with initial stress and fibre
    directions.

    The law is history-independent and does not store any internal state,
    but it is fully compatible with the mechanicalConstitutiveLaw framework
    and may be used alongside history-dependent laws in multi-material
    simulations.

    Material parameters are provided in terms of Young's modulus and
    Poisson's ratio, specified as dimensioned quantities and validated at
    construction time.

    The implementation is intended as a simple, robust finite-strain
    reference law and as a demonstrator for the finite-strain constitutive
    infrastructure within solids4foam.

SourceFiles
    MooneyRivlinElasticMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef MooneyRivlinElasticMechanicalConstitutiveLaw_H
#define MooneyRivlinElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class MooneyRivlinElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class MooneyRivlinElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private Data

        //- Density
        const dimensionedScalar rho_;

        //- Young's modulus
        dimensionedScalar c10_;

        //- Mooney-Rivlin coefficient
        dimensionedScalar c01_;

        //- Mooney-Rivlin coefficient
        dimensionedScalar c11_;

        //- Young's modulus, only used to derive K from nu
        dimensionedScalar E_;

        //- Poisson's ration
        dimensionedScalar nu_;

        //- First Lamé parameter
        dimensionedScalar lambda_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;


    // Private Member Functions

        //- No copy construct
        MooneyRivlinElasticMechanicalConstitutiveLaw
        (
            const MooneyRivlinElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const MooneyRivlinElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("MooneyRivlinElastic");

    // Constructors

        //- Construct from dictionary
        explicit MooneyRivlinElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~MooneyRivlinElasticMechanicalConstitutiveLaw() = default;


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
            return kappa_;
        }

    // Edit

        //- Small-strain constitutive evaluation
        virtual void evaluate
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
