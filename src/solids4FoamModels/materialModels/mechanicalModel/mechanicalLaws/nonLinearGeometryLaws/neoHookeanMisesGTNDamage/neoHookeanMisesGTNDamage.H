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
    neoHookeanMisesGTNDamage

Description
    Large strain elasto-plastic model with Gurson-Tveergard-Needleman (GTN) micro-mechanical damage model. This class also provides users the option to include shear effects and the non-local porosity f.

    Large strain elasto-plastic model with GTN yield equation uses algorithmic approach similar to those outined in Achouri et Al. (Doi: 10.1016/j.commatsci.2013.01.035) and Leclerc et Al. (Doi: 10.1016/j.ijplas.2019.11.010).

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    The plastic hardening law is defined in a user file (see plasticStrainVsYieldStress in tutorial case) 

    For GTN parameters, the user can specify:
        - the yield equation paramters q1 q2 and q3
        - void nucleation paramters epsilonN, sN, fN
        - void coalescence parameters f0, fC, fU and fF;
        - shear parameter kW
    
    For non-local porosity the user can speicfy 

SourceFiles
    neoHookeanMisesGTNDamage.C

Author
    Andrew Whelan, UCD.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMisesGTNDamage_H
#define neoHookeanMisesGTNDamage_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class neoHookeanMisesGTNDamage Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMisesGTNDamage
:
    public mechanicalLaw
{
    // Private data

        const scalar relaxationFactor_;

        //Switch to include shear effects in the evolution of the porosity
        Switch includeShear_;

        //Switch to include coalescene of voids effects
        Switch includeCoalescence_;

        //characteristic length for non-local porosity
        dimensionedScalar charLength_;

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

        //GTN yield equaiton parameters
        dimensionedScalar q1_;
        dimensionedScalar q2_;
        dimensionedScalar q3_;

        //initial porosity
        dimensionedScalar f0_;

        //parameters for the nucleation of voids
        dimensionedScalar fN_;
        dimensionedScalar epsilonN_;
        dimensionedScalar sN_;

        //coalescence of voids parameters
        dimensionedScalar fC_;
        dimensionedScalar fU_;
        dimensionedScalar fF_;

        //shear parameter
        dimensionedScalar kw_;

        //- Relative deformation gradient volField
        volTensorField* relFPtr_;

        //- Jacobian of the toal deformation gradient volField
        volScalarField* JPtr_;

        //- Total deformation gradient volField
        volTensorField* FPtr_;

        //- Table of yield (Cauchy) stress versus true plastic strain
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //porosity
        volScalarField f_;

        //effective porosity
        volScalarField fStar_;

        //non-local porosity
        volScalarField fNonLocal_;

        //gradient of the non-local porosity
        volVectorField gradfNonLocal_;

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
        neoHookeanMisesGTNDamage(const neoHookeanMisesGTNDamage&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanMisesGTNDamage&);

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
        volTensorField& F1();

        //calculate porosity
        void calculatef
        (
            const scalar& DHyd,
            const scalar& DqEpsilonP,       
            scalar& f,  
            scalar& fStar,          
            const scalar& fOld, 
            const scalar& fNonLocal,               
            const scalar& epsilonPEqOld,    
            const scalar& DEpsilonPEq,
            const symmTensor& DEpsilonP,
            const symmTensor& devT,
            const scalar& p
        );

        //- Iteratively calculate deviatoric, spherical and equivalent plastic strain
        void smallStrainReturnMap
        (
            const scalar& pTrial, 
            const symmTensor& sTrial, 
            const scalar& qTrial,                         
            scalar& DHyd,
            scalar& DqEpsilonP,       
            const scalar& fNonLocal,             
            const scalar& fStar,               
            const scalar& epsilonPEqOld,    
            scalar& DEpsilonPEq,
            symmTensor& plasticN
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
    TypeName("neoHookeanMisesGTNDamage");

    // Static data members


    // Constructors

        //- Construct from dictionary
        neoHookeanMisesGTNDamage
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~neoHookeanMisesGTNDamage();


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
