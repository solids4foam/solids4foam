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
    coupledPressureDisplacementSolid

Description
    Block-coupled pressure-displacement solid model for incompressible and
    nearly incompressible solids.

    The model solves displacement increment and hydrostatic pressure increment
    in a coupled system and is currently available for foam-extend builds.

Author
    Zeljko Tukovic, FSB Zagreb.

SourceFiles
    coupledPressureDisplacementSolid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPressureDisplacementSolid_H
#define coupledPressureDisplacementSolid_H

#include "solidModel.H"

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class coupledPressureDisplacementSolid Declaration
\*---------------------------------------------------------------------------*/

class coupledPressureDisplacementSolid
:
    public solidModel
{
    // Private data

        //- Switch to enable/disable nonlinear geometry
        Switch nonLinear_;

        //- Block vector field for velocity and pressure
        volVector4Field DDDp_;

        //- Velocity increment field
        volVectorField DU_;

        //- Hydrostatic pressure increment field
        volScalarField Dp_;

        //- Hydrostatic pressure field
        volScalarField p_;

        //- Hydrostatic pressure increment field
        surfaceScalarField Dpf_;

        //- Hydrostatic pressure field
        surfaceScalarField pf_;

        //- Hydrostatic pressure field gradient
        volVectorField gradDp_;

        //- Surface displacement field
        surfaceVectorField Df_;

        //- Surface velocity field
        surfaceVectorField Uf_;

        //- Surface velocity increment field
        surfaceVectorField DUf_;

        //- Surface displacement field
        surfaceVectorField DDf_;

        //- Surface area field
        surfaceVectorField curSf_;

        //- Diagonal coefficient from laplacian term
        volScalarField AD_;

        //- Diagonal coefficient from temporal term
        volScalarField ADt_;

        //- Surface flux
        surfaceScalarField phi_;

        //- Non-linear force correction
        surfaceVectorField nonLinForceCorrection_;

        //- Non-linear force
        surfaceVectorField force_;

        //- Solid-interface traction for postprocessing
        volVectorField solidInterfaceTraction_;

        //- Acceleration
        volVectorField A_;

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Surface displacement gradient
        surfaceTensorField gradDDf_;

        //- Surface displacement gradient
        surfaceTensorField gradDf_;

        //- Surface displacement gradient
        surfaceVectorField tGradDDn_;

        //- Surface displacement gradient
        surfaceVectorField nGradDDn_;

        //- Implicit stiffness surface field
        mutable surfaceScalarField impKf_;

        //- Inverse bulk modulus
        volScalarField rKappa_;

        //- lambda/kappa = 3*nu/(1+nu)
        volScalarField nuCoeff_;

        //- Debug switch for the print-out of residuals
        const Switch debug_;

        //- Viscosity/damping coefficient
        mutable dimensionedScalar K_;

        //- Relative tolerance
        const scalar relativeTol_;

        //- Reference configuration mesh points;
        const pointIOField refPoints_;

        //- Reference configuration mesh points;
        // TODO: move to separate class
        Switch consistentRhieChow_;

        //- Are we in the reference configuration?
        bool refConfig_;

        //- Composite temporal scheme (temporary)
        Switch composite_;

        //- Composite temporal scheme;
        // Switch skewCorr_;

        //- Use standard procedure for displ grad calculation
        Switch stdDispGrad_;

        //- Temporal discretisation coefficient;
        dimensionedScalar coeff_;

        //- Owner to neighbour cell centre connection vector
        surfaceVectorField delta_;

        //- Weighting factor
        surfaceScalarField weights_;

        //- Convergence data file
        autoPtr<OFstream> convergDataFilePtr_;

        //- Switch to write the residual field
        const Switch writeConvergData_;

        //- Surface displacement increment field
        mutable surfaceVectorField* interpDDfPtr_;

        //- Surface velocity increment field
        mutable surfaceVectorField* interpDUfPtr_;

        //- Surface displacement field
        mutable surfaceVectorField* interpDfPtr_;

        //- Surface velocity field
        mutable surfaceVectorField* interpUfPtr_;


    // Private Member Functions

        surfaceVectorField& interpDDf()
        {
            return *interpDDfPtr_;
        }

        surfaceVectorField& interpDUf()
        {
            return *interpDUfPtr_;
        }

        surfaceVectorField& interpDf()
        {
            return *interpDfPtr_;
        }

        surfaceVectorField& interpUf()
        {
            return *interpUfPtr_;
        }

        //- Return the displacement residual
        scalar residual(const volVectorField& D) const;

        //- Disallow default bitwise copy construct
        coupledPressureDisplacementSolid
        (
            const coupledPressureDisplacementSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const coupledPressureDisplacementSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            if (nonLinear_)
            {
                return nonLinearGeometry::TOTAL_LAGRANGIAN;
            }
            else
            {
                return nonLinearGeometry::LINEAR_GEOMETRY;
            }
        }


public:

    //- Runtime type information
    TypeName("coupledPressureDisplacementSolid");

    // Constructors

        //- Construct from components
        coupledPressureDisplacementSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~coupledPressureDisplacementSolid()
        {}

    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for DD
                return DD();
            }

            //- This solid model solves for DD instead of D
            virtual bool incremental() const
            {
                return true;
            }

            //- Return nonLinear geometry switch
            const Switch& nonLinear() const
            {
                return nonLinear_;
            }

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Traction boundary hydrostatic pressure
            tmp<scalarField> tractionBoundaryHydPressure
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Move mesh to current configuration
            void moveToCurConfig(const bool prevConfig = false);

            //- Move mesh to reference configuration
            void moveToRefConfig();

            //- Set traction at specified patch
            virtual void setTraction
            (
                fvPatchVectorField& tractionPatch,
                const vectorField& traction
            );

            //- Set pressure at specified patch
            virtual void setPressure
            (
                fvPatchVectorField& pressurePatch,
                const scalarField& pressure
            );

            virtual void setTraction
            (
                const label interfaceI,
                const label patchID,
                const vectorField& faceZoneTraction
            );

            //- Face zone acceleration
            virtual tmp<vectorField> faceZoneAcceleration
            (
                const label interfaceI
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
