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
    IQNILSCouplingInterface

Description
    Strong IQN-ILS Dirichlet-Neumann coupling, as developed by Degroote et al.

    A fluid solid interaction solver with IQN-ILS coupling algorithm.
    J. Degroote, K.-J. Bathe and J. Vierendeels.
    Performance of a new partitioned procedure versus a monolithic
    procedure in fluid-solid interaction. Computers & Solids

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    IQNILSCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef IQNILSCouplingInterface_H
#define IQNILSCouplingInterface_H

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class IQNILSCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class IQNILSCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- Fixed under-relaxation factor
        const scalar relaxationFactor_;

        //- Coupling reuse
        label couplingReuse_;

        //- Absolute tolerance for filtering repeated/near-dependent modes
        const scalar minSignificant_;

        //- Relative tolerance for QR2-style filtering of near-dependent modes
        const scalar relMinSignificant_;

        //- Normalize IQN-ILS secant columns before the QR solve
        const bool normalizeCouplingColumns_;

        //- Relative tolerance used to drop nearly singular directions in the
        //  QR back-substitution
        const scalar qrSolveTolerance_;

        //- Relative QR-stage filter tolerance used to ignore weak older secant
        //  columns for the current solve without deleting history
        const scalar qrFilterTolerance_;

        //- Build one combined IQN-ILS least-squares system across all
        //  interfaces, with columns aligned by recency
        const bool combinedCouplingSystem_;

        //- Build the IQN-ILS least-squares system with a preCICE-like QR2 pass
        const bool preciceStyleCouplingQR_;

        //- Keep previous time-step modes for the first outer iteration even
        //  when couplingReuse is zero, then drop them before iteration two
        const bool reusePreviousStepModes_;

        //- Minimum number of cached secant modes required before previous-step
        //  reuse is considered ready
        const label minPreviousStepModesForReuse_;

        //- Limit the norm of the first reused IQN-ILS correction relative to
        //  the current interface residual norm. Negative disables the limit.
        const scalar maxReuseUpdateNormRatio_;

        //- Perform a second MGS pass while building the QR factorization
        const bool reorthogonalizeCouplingColumns_;

        //- Enable residual-sum preconditioning
        const bool residualSumPreconditioning_;

        //- Small regularization added to the inverse residual magnitude
        const scalar residualPreconditioningEpsilon_;

        //- Predict solid
        const bool predictSolid_;

        //- List of coupling field V
        List<DynamicList<vectorField> > fluidPatchesPointsV_;

        //- List of coupling field W
        List<DynamicList<vectorField> > fluidPatchesPointsW_;

        //- List of secondary coupling field W for mapped solid tractions
        List<DynamicList<vectorField> > solidPatchesFacesTractionW_;

        //- List of coupling field T
        List<DynamicList<scalar> > fluidPatchesPointsT_;

        //- Cached secant information from the previous converged time step
        List<DynamicList<vectorField> > previousStepPatchesPointsV_;

        //- Cached secant updates from the previous converged time step
        List<DynamicList<vectorField> > previousStepPatchesPointsW_;

        //- Cached secondary traction secants from the previous converged step
        List<DynamicList<vectorField> > previousStepPatchesFacesTractionW_;

        //- Backup secants from the latest converged time step that contained
        //  usable IQN-ILS modes. These are used if the immediately previous
        //  time step converged in one iteration and produced no new modes.
        List<DynamicList<vectorField> > previousNonEmptyStepPatchesPointsV_;

        //- Backup secant updates from the latest converged time step that
        //  contained usable IQN-ILS modes
        List<DynamicList<vectorField> > previousNonEmptyStepPatchesPointsW_;

        //- Backup secondary traction secants from the latest converged time
        //  step that contained usable IQN-ILS modes
        List<DynamicList<vectorField> >
            previousNonEmptyStepPatchesFacesTractionW_;


    // Private Member Functions

        //- Reuse coupling
        label couplingReuse() const;

        //- Remove a stored coupling mode for the given interface
        void removeCouplingMode(const label interfaceI, const label modeI);

        //- Remove modes older than the specified minimum time index
        void removeOldCouplingModes
        (
            const label interfaceI,
            const label minTimeIndexToKeep
        );

        //- Clear all active coupling modes for the given interface
        void clearCouplingModes(const label interfaceI);

        //- Cache the current time-step coupling modes for next-step reuse
        void cacheCurrentStepModes();

        //- True if the previous time-step cache should be used in the current
        //  IQN-ILS solve for this interface
        bool usePreviousStepModesInSolve(const label interfaceI) const;

        //- True if the latest non-empty cached time step should be used
        //  because the immediately previous time step converged without
        //  generating secant modes
        bool usePreviousNonEmptyStepModesInSolve(const label interfaceI) const;

        //- Number of secant modes available to the current IQN-ILS solve,
        //  combining active and optionally cached previous-step modes
        label nSolveCouplingModes(const label interfaceI) const;

        //- Return the nth newest secant mode available to the current
        //  IQN-ILS solve. Active modes are preferred over previous-step modes.
        void solveCouplingMode
        (
            const label interfaceI,
            const label newestModeI,
            vectorField& v,
            vectorField& w
        ) const;

        //- Return the nth newest secondary traction secant update available to
        //  the current IQN-ILS solve.
        void solveSecondaryTractionMode
        (
            const label interfaceI,
            const label newestModeI,
            vectorField& wTraction
        ) const;

        //- Total number of interface points across all interfaces
        label combinedCouplingSystemSize();

        //- Number of aligned secant columns in the combined IQN-ILS system
        label nCombinedCouplingModes();

        //- Return the nth newest aligned secant column in the combined system
        void combinedCouplingMode
        (
            const label newestModeI,
            vectorField& v,
            vectorField& w
        );

        //- Return the combined interface residual for all interfaces
        vectorField combinedResidual();

        //- Set all interface displacements from a combined correction vector
        void applyCombinedCorrection(const vectorField& correction);

        //- Filter repeated/near-dependent coupling modes
        void filterCouplingModes(const label interfaceI);

        //- True if the preCICE-style QR2 solve should be used for this
        //  interface in the current outer iteration
        bool usePreciceStyleCouplingQR(const label interfaceI) const;

        //- True if the combined LS system should use the preCICE-style QR2
        //  solve in the current outer iteration
        bool useCombinedPreciceStyleCouplingQR() const;

        //- Filter near-dependent modes using a relative QR2-style criterion
        void filterRelativeCouplingModes(const label interfaceI);

        //- Orthogonalize and normalize a new QR column using a preCICE-like
        //  repeated MGS procedure. Returns -1 if the column should be discarded.
        label orthogonalizeCouplingColumn
        (
            vectorField& v,
            scalarField& r,
            scalar& rho,
            const DynamicList<vectorField>& qBasis
        ) const;

        //- Build the IQN-ILS QR system using newest-first QR2-style insertion
        label buildPreciceStyleCouplingQR
        (
            const label interfaceI,
            DynamicList<vectorField>& solveQ,
            DynamicList<vectorField>& solveW,
            RectangularMatrix<scalar>& R
        ) const;

        //- Build a preCICE-style QR2-filtered combined IQN-ILS system
        label buildPreciceStyleCombinedCouplingQR
        (
            DynamicList<vectorField>& solveQ,
            DynamicList<vectorField>& solveW,
            RectangularMatrix<scalar>& R
        );

        //- Disallow default bitwise copy construct
        IQNILSCouplingInterface(const IQNILSCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const IQNILSCouplingInterface&);


public:

    //- Runtime type information
    TypeName("IQNILS");


    // Constructors

        //- Construct from components
        IQNILSCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~IQNILSCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Calculate interface displacement
            virtual void updateDisplacement();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
