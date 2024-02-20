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
    neoHookeanMisesPhaseFieldDamage

Description
    This law is ported from https://github.com/awhelanUCD/solids4foam where it
    is called logStrainPhaseFieldFracture.

    A ductile phase-field fracture model coupled with elasto-plastic model. This
    class gives user the option to incorporate triaxiality effects into the
    phase-field model. The phase field is described by the variable "d" with
    initial (unbroken/undamaged) value of 1.0.

    Chacteristics
        - Neo-Hookean elasticity with Mises/J2 plasticity, Etrovic & Bathe,
          DOI:10.1002/nme.162030060.
        - Ductile phase-field fracure model, Borden et al.,
          DOI:10.1016/j.cma.2016.09.005)

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu), or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a XY interpolate file, e.g.
    see "plasticStrainVsYieldStress" in tutorial case.

    For the phase-field fracture model, the user can specify:
        - Critical fracture energy (Gc), characteristic length (lc), plastic
          work threshold threshold (w0)
        - Triaxiality Parameters (d1, d2, d3)

SourceFiles
    neoHookeanMisesPhaseFieldDamage.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMisesPhaseFieldDamage_H
#define neoHookeanMisesPhaseFieldDamage_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class neoHookeanMisesPhaseFieldDamage Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMisesPhaseFieldDamage
:
    public mechanicalLaw
{
    // Private data

        //- Characteristic length for non-local gradient damage equation
        dimensionedScalar charLength_;

        //- Critical fracture energy
        dimensionedScalar Gc_;

        //- Triaxiality parameters
        dimensionedScalar d1_;
        dimensionedScalar d2_;
        dimensionedScalar d3_;

        //- Parameters for weighting elastic and plastic energy contributions,
        //  generally not used
        const scalar betaE_;
        const scalar betaP_;

        //- Plastic work threshold
        dimensionedScalar w0_;

        //- Solve the phase field equation monolithically, as opposed to the
        //  staggered approach
        Switch monolithic_;

        //- Incorporate triaxiality effects into the phase field fracture field
        Switch triaxiality_;

        //- Initial density
        const dimensionedScalar rho_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Relative deformation gradient volField
        volTensorField* relFPtr_;

        //- Jacobian of the toal deformation gradient volField
        volScalarField* JPtr_;

        //- Total deformation gradient volField
        volTensorField* FPtr_;

        //- Table of post-yield stress versus plastic strain
        //  For large strain total Lagrangian models, it expects
        //  the true stress and true plastic strain.
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Phase field variable d
        volScalarField Dpf_;

        //- Dimensionless unity field
        const volScalarField one_;

        //- Gradient of the phase field variable
        volVectorField gradDpf_;

        //- Elastic strain energy
        volScalarField strainEnergy_;

        //- History variable for the elastic strain energy
        volScalarField HStrainEnergy_;

        //- Plastic strain energy
        volScalarField plasticStrainEnergy_;

        //- Deviatoric elastic strain energy
        volScalarField devStrainEnergy_;

        //- Volumetric elastic strain energy
        volScalarField volStrainEnergy_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        //- Switch to enable solution of a pressure equation, which ensures a
        //  smooth sigmaHyd field
        Switch smoothPressure_;

        //Kirchhoff stress tensor
        volSymmTensorField tau_;

        //- Cauchy yield stress
        volScalarField sigmaY_;

        //- Incremental change in sigmaY
        volScalarField DSigmaY_;

        //elastic strain
        volSymmTensorField epsilonEl_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        volSymmTensorField bEbarTrial_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        volSymmTensorField bEbar_;

        //- Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;

        //- Equivalent plastic strain
        volScalarField epsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        //- An iterative procedure is used when the yield stress is a nonlinear
        //  function of plastic strain
       // const bool nonLinearPlasticity_;

        //- Update bEbar consistently with the assumption that det(bEbar) == 1
        //  defaults to off
        const Switch updateBEbarConsistent_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        neoHookeanMisesPhaseFieldDamage(const neoHookeanMisesPhaseFieldDamage&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanMisesPhaseFieldDamage&);

        symmTensor logm(const symmTensor& T);

        //- Exponential of a symmTensor
        symmTensor expm(const symmTensor& T);

        //- Make the relF field
        void makeRelF();

        //- Return a reference to the relF field
        volTensorField& relF();

        //- Make the J field
        void makeJ();

        //- Return a reference to the J field
        volScalarField& J();

        //- Make the F field
        void makeF();

        //- Return a reference to the F field
        volTensorField& F();

        //- Iteratively calculate plastic multiplier increment (DLambda)
        //  and current yield stress using Newton's method

        void smallStrainReturnMap
        (
            const symmTensor& Ee,
            scalar& sigmaY,
            scalar& DSigmaY,
            scalar& DLambda,
            const scalar& epsilonPEqOld,
            scalar& activeYield,
            const scalar& gc
        ) const;

        //- Calcualte Ibar such that det(bEbar) == 1
        tmp<volScalarField> Ibar
        (
            const volSymmTensorField& devBEbar
        );

        //- Calcualte Ibar such that det(bEbar) == 1
        tmp<surfaceScalarField> Ibar
        (
            const surfaceSymmTensorField& devBEbar
        );

public:

    //- Runtime type information
    TypeName("neoHookeanMisesPhaseFieldDamage");

    // Static data members


    // Constructors

        //- Construct from dictionary
        neoHookeanMisesPhaseFieldDamage
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~neoHookeanMisesPhaseFieldDamage();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        virtual void calcPhase();

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();

        //- Return the desired new time-step
        virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
