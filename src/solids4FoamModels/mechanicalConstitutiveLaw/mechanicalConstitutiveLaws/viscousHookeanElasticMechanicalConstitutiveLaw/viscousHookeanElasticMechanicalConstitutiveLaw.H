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
    Foam::viscousHookeanElasticMechanicalConstitutiveLaw

Description
    Linear viscoelastic constitutive law: Hookean elasticity in the
    volumetric response with a generalised Maxwell (Prony series) deviatoric
    response.

    The deviatoric stress relaxes through n Maxwell arms,

        s   = 2 mu dev(symm(gradD))
        h_i = exp(-dt/tau_i) h_i,0 + exp(-dt/(2 tau_i)) (s - s_0)
        sigma = k tr(gradD) I + gammaInf s + sum_i gamma_i h_i

    following Simo & Hughes (1998), eq. 10.3.12. The instantaneous modulus is
    E0 = EInfinity + sum(E_i), with gammaInf = EInfinity/E0 and
    gamma_i = E_i/E0.

    This is the first law to take a value from the inputs object: the update
    is a function of the time increment, and so is the scalar tangent, whose
    relaxation factor is gammaInf + sum_i gamma_i exp(-dt/(2 tau_i)). The time
    increment is not a deformation measure and no longer lives in the
    kinematics.

    The Williams-Landel-Ferry temperature shift of the relaxation times is
    supported in the algorithm and asks for the temperature as a live input.
    Nothing supplies that input yet, so a case enabling it fails with a clear
    message rather than silently ignoring the shift; it will work unchanged
    once the manager can scatter a field into the inputs.

    Elastic parameters are specified as EInfinity, nu, a list of moduli E and
    a matching list of relaxationTimes.

    References:
      - J.C. Simo & T.J.R. Hughes, Computational Inelasticity, Springer, 1998.

SourceFiles
    viscousHookeanElasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef viscousHookeanElasticMechanicalConstitutiveLaw_H
#define viscousHookeanElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "scalarList.H"
#include "Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class viscousHookeanElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class viscousHookeanElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private Data

        //- Density
        const dimensionedScalar rho_;

        //- Long-term (equilibrium) Young's modulus
        const dimensionedScalar EInf_;

        //- Moduli of the Maxwell arms
        const scalarList E_;

        //- Relaxation times of the Maxwell arms
        const scalarList tau_;

        //- Poisson's ratio
        const dimensionedScalar nu_;

        //- Relative equilibrium modulus
        dimensionedScalar gammaInf_;

        //- Relative moduli of the Maxwell arms
        scalarList gamma_;

        //- Shift the relaxation times with temperature
        const Switch williamsLandelFerryShift_;

        //- Williams-Landel-Ferry coefficients
        const dimensionedScalar C1_;
        const dimensionedScalar C2_;
        const dimensionedScalar Tref_;

        //- First Lamé parameter
        dimensionedScalar lambda_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;

        //- State field name of the i-th Maxwell arm
        static word arm(const label i);


    // Private Member Functions

        //- No copy construct
        viscousHookeanElasticMechanicalConstitutiveLaw
        (
            const viscousHookeanElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const viscousHookeanElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("viscousHookeanElastic");

    // Constructors

        //- Construct from dictionary
        explicit viscousHookeanElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~viscousHookeanElasticMechanicalConstitutiveLaw() = default;


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
        //- Initialise the relaxing stress state
        virtual void initialiseState
        (
            mechanicalConstitutiveLawState& state
        ) const;

        virtual void evaluate
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
