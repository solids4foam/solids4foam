/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    neoHookeanMisesLemaitreDamage

Description
    This law is ported from https://github.com/awhelanUCD/solids4foam where it
    is called logStrainLemaitreDamage.

    Lemaitre damage model coupled with elasto-plastic model. This class also give
    users the option to incorporate non-local gradient damage and to use the
    Malcher variant of the Lemaitre damage law.

    Chacteristics
     - Neo-Hookean elasticity with Mises/J2 plasticity, Etrovic & Bathe,
       DOI: 10.1002/nme.162030060.
     - Lemaitre damage model with crack closure effects, Lemaitre,
       ISBN: 978-3-540- 21503-5.
     - non local gradient damage, Peerlings et al.,
       DOI: 10.1016/S0020-7683(01)00087-7

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu), or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a XY interpolate file, e.g.
    see "plasticStrainVsYieldStress" in tutorial case. 

    For Lemaitre damage parameters, the user can specify:
        - Damage denominator (s0)
        - Damage exponent (b),
        - Plastic strain threshold (epsilonD)

    For non-local gradient damage model, the user can specify:
        - The characteristic length (lc) (for local damage model lc=0)

    For the Malcher variant the user can specify:
        - Damage denominators s00 and s33

  
SourceFiles
    neoHookeanMisesLemaitreDamage.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD. 
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMisesLemaitreDamage_H
#define neoHookeanMisesLemaitreDamage_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
#include "Eigen/Dense"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class neoHookeanMisesLemaitreDamage Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMisesLemaitreDamage
:
    public mechanicalLaw
{
    // Private data

        const scalar relaxationFactor_;

        const string lemaitreType_;

        //- Table of post-yield stress versus plastic strain
        //  For large strain total Lagrangian models, it expects
        //  the true stress and true plastic strain.
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Characteristic length for non-local damage
        dimensionedScalar charLength_;

        //- Plastic strain threshold
        dimensionedScalar epsilonD_;

        //- Initial density
        const dimensionedScalar rho_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;
      
        //- Lemaitre exponent
        dimensionedScalar b_;

        //- Lemaitre denominator
        dimensionedScalar s0_;

        //- Lemaitre denominator for Malcher variant
        dimensionedScalar s00_;

        //- Lemaitre denominator for Malcher variant
        dimensionedScalar s33_;

        //- Critical damage parameter
        scalar damageC_;

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

        //damage
        volScalarField damage_;

        //non-local damage
        volScalarField damageNonLocal_;

        //gradient of the non-local damage
        volVectorField gradDamageNonLocal_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigma)/3.0 == -p,
        //  where p is the hydrostatic pressure
        volScalarField sigmaHyd_;

        //- Switch to enable solution of a pressure equation, which ensures a
        //  smooth sigmaHyd field
        Switch smoothPressure_;

        //Kirchhoff stress
        volSymmTensorField tau_;

        //- Cauchy yield stress
        volScalarField sigmaY_;

        //- Incremental change in sigmaY
        volScalarField DSigmaY_;

        volSymmTensorField epsilonEl_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        volSymmTensorField sigmaPositive_;

        volSymmTensorField sigmaNegative_;

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

        volScalarField triaxiality_;

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
        neoHookeanMisesLemaitreDamage(const neoHookeanMisesLemaitreDamage&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanMisesLemaitreDamage&);

        //log of a symmTensor
        symmTensor logm(const symmTensor& T);

        //- Exponential of a symmTensor
        symmTensor expm(const symmTensor& T);


        //decompose SymmTensor into it's postive and negative components
        void decompose
        (
            const symmTensor& sigma,
            symmTensor& sigmaPositive,
            symmTensor& sigmaNegative
        );

        //macauley bracket
        scalar macauley
        (
            scalar x
        );

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

        //calculate the normalized third invariant of the deviatoric stress
        scalar calcZheta(const symmTensor& S);

        // calculate the damage
        void damage
        (
            const scalar& DLambda,
            const scalar& DEpsilonPEq,              
            const scalar& epsilonPEqOld,  
            const scalar& triaxiality,                  
            scalar& damage, 
            scalar& damageOld,
            const scalar&  damageNonLocal,            
            const symmTensor& tau
        );

        // small strain return-map to calculate plastic increment
        void smallStrainReturnMap
        (
            const symmTensor& Ee,  
            scalar& sigmaY,                      
            scalar& DLambda, 
            scalar& damage, 
            const scalar& damageOld,       
            const scalar& damageNonLocal,                           
            const scalar& epsilonPEqOld,
            scalar& activeYield
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
    TypeName("neoHookeanMisesLemaitreDamage");

    // Static data members


    // Constructors

        //- Construct from dictionary
        neoHookeanMisesLemaitreDamage
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~neoHookeanMisesLemaitreDamage();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

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
