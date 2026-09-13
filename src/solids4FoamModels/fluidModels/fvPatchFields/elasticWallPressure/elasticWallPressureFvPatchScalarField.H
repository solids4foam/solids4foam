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
    elasticWallPressureFvPatchScalarField

Description
    Robin boundary condition for pressure of a flow in a fluid-solid
    interaction problem which accounts for the solid inertia.

    The Robin coefficient rho_s*hs is an estimate of the solid interface
    impedance (added mass per unit area) and controls the convergence of the
    Robin-Neumann coupling: the iteration contracts fastest when rho_s*hs
    matches the solid impedance, slows when it is smaller, and diverges for
    strong added-mass problems when it is more than about twice as large.

    The solid virtual thickness hs is selected with hsModel:
    - pWaveSpeed: hs = c_p*dt_eff, the impedance of an elastic half-space,
      where c_p = sqrt(impK/rho_s) (the previous default);
    - constant: hs = constantHs, or a nonuniform "hs" entry on the patch;
    - thicknessLimited: hs = l*tanh(H/l) (or min(l, H)), with l = c_p*dt_eff
      and H the local solid wall thickness found by ray casting through the
      solid; H is halved when the ray exits through another FSI interface,
      i.e. for structures wetted on both sides;
    - secant (default): starts from seedModel (thicknessLimited) and rescales
      the coefficient using the solid and fluid interface impedances measured
      from the changes in interface pressure and solid acceleration between
      FSI iterations. By default the scale is updated every iteration
      (secantUpdate iteration) to minimise the largest predicted contraction
      factor over the recently observed modes (secantFit minimax).

    thicknessLimited is a safe a-priori alternative (a lower bound of the
    solid impedance), while pWaveSpeed can diverge for thin or nearly
    incompressible walls and large time steps. With writeDiagnostics, the
    measured impedances and the predicted and observed contraction factors
    are written to postProcessing/robinCoefficient_<patch>.dat. See
    tutorials/fluidSolidInteraction/robinCouplingStudy for the study behind
    these choices.

    Coupling: the Robin-Neumann coupling is designed for unrelaxed fixed-point
    FSI iterations, i.e. fluidSolidInterface fixedRelaxation with
    relaxationFactor 1, which is the recommended setting. With Aitken, IQN-ILS
    or a relaxation factor below 1 the fluid mesh lags the solid within the
    iterations: a warning is issued, the interface flux is not made consistent
    with the mesh motion and the leakage is only reported.

    Kinematic consistency: the Robin condition uses the solid acceleration of
    the previous FSI iteration, as in Tukovic et al. (2018). At convergence
    the interface flux of that formulation equals the mesh (swept-volume)
    flux only in the continuum; pimpleFluid and interFluid (OpenFOAM.com)
    therefore build the explicit part of the interface flux from the mesh
    flux (robinKinematicConsistency) so that the converged flux equals it
    exactly, as for a Dirichlet velocity condition, without changing the
    Robin condition itself.

    Example:
    \verbatim
    interface
    {
        type            elasticWallPressure;
        hsModel         secant;             // Optional, default
        hsScale         1;                  // Optional
        value           uniform 0;
    }
    \endverbatim

    Reference:
    Tuković et al. Added Mass Partitioned Fluid–Structure Interaction Solver
    Based on a Robin Boundary Condition for Pressure. OpenFOAM. 2019;1–22.

SourceFiles
    elasticWallPressureFvPatchScalarField.C
    elasticWallPressureFvPatchScalarFieldRobinCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef elasticWallPressureFvPatchScalarField_H
#define elasticWallPressureFvPatchScalarField_H

#include "robinFvPatchFields.H"
#include "Switch.H"
#include "OFstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidSolidInterface;

/*---------------------------------------------------------------------------*\
                 Class elasticWallPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class elasticWallPressureFvPatchScalarField
:
    public robinFvPatchField<scalar>
{
    // Private data

        //- Pressure applied at the solid side of the interface
        //  in the previous fsi iteration
        scalarField prevPressure_;

        //- Acceleration of the solid side of the interface
        //  in the previous fsi iteration
        vectorField prevAcceleration_;

        //- A pointer to the solid density times the solid virtual thickness
        //  field mapped to the fluid patch, i.e. the Robin coefficient
        mutable autoPtr<scalarField> rhoSolidHsPtr_;

        //- Does the fluid model build the interface flux such that the
        //  converged flux equals the mesh flux (no leakage)
        bool kinematicConsistency_;

        //- Is the fluid interface moved to the solid interface in every FSI
        //  iteration (set by the FSI coupling)
        bool fluidMeshFollowsSolid_;

        //- Time index at which the Robin number was last reported
        label robinNumberTimeIndex_;

        //- Optional constant solid thickness
        //  If constantHs < SMALL, it indicates that it should not be used and
        // instead a virtual solid thickness is calculated on the fly
        const scalar constantHs_;

        //- Settings for the Robin coefficient (rho_s*hs) estimate
        struct hsControls
        {
            //- pWaveSpeed, constant, thicknessLimited or secant
            word model;

            //- Was hsModel given explicitly
            bool modelSpecified;

            //- Scale factor applied to the estimated coefficient
            scalar scale;

            //- Blending for thicknessLimited: tanh or min
            word thicknessBlend;

            //- Halve the thickness for walls wetted on both sides
            Switch twoSidedHalving;

            //- Optional wave speed override; the impK wave speed is used if
            //  negative
            scalar waveSpeed;

            //- A-priori model used to seed the secant model
            word seedModel;

            //- Secant update frequency: timeStep or iteration
            word secantUpdate;

            //- Secant fit: median (of the per-iteration solid impedance
            //  fits) or minimax (of the predicted contraction factor over
            //  the recorded solid/fluid impedance samples)
            word secantFit;

            //- Fit a spatially varying secant scale factor (experimental; the
            //  smoothing only uses processor-local face neighbours)
            Switch secantSpatial;

            //- Number of patch smoothing sweeps for the spatial secant fit
            label secantSmoothing;

            //- Bounds on the secant scale factor relative to the seed
            scalar secantMinFactor;
            scalar secantMaxFactor;

            //- Maximum decrease factor of the secant scale per update
            scalar secantMaxChange;

            //- Maximum increase factor of the secant scale per update;
            //  increases are limited more than decreases because an
            //  overestimated coefficient can diverge
            scalar secantMaxIncrease;

            //- Number of previous time steps whose impedance samples are
            //  kept for the minimax fit
            label secantMemory;

            //- Halve the coefficient when the FSI iterations diverge
            Switch divergenceSafeguard;

            //- Write per-iteration diagnostics
            Switch writeDiagnostics;

            //- Construct with default settings
            hsControls();

            //- Construct from dictionary
            hsControls(const dictionary& dict, const scalar constantHs);

            //- Write the non-default settings
            void write(Ostream& os) const;
        };

        //- Robin coefficient settings
        hsControls hsControls_;

        //- Optional user-specified nonuniform hs on the fluid patch
        autoPtr<scalarField> hsUserPtr_;

        //- Local solid wall thickness on the solid patch
        mutable autoPtr<scalarField> solidThicknessPtr_;

        //- A-priori (seed) coefficient rho_s*hs on the fluid patch, before
        //  hsScale is applied
        mutable autoPtr<scalarField> seedCoeffPtr_;

        //- Time index at which seedCoeffPtr_ was last updated
        mutable label seedCoeffTimeIndex_;

        //- Secant scale factor field applied to the seed coefficient
        scalarField secantScale_;

        //- Scale factor from the divergence safeguard
        scalar safeguardScale_;

        //- Time index of the secant history
        label secantTimeIndex_;

        //- Number of FSI iterations recorded in the current time step
        label secantIter_;

        //- Interface pressure at the previous FSI iteration
        scalarField iterPressure_;

        //- Solid normal acceleration at the previous FSI iteration
        scalarField iterSolidAccN_;

        //- Fluid normal acceleration at the previous FSI iteration
        scalarField iterFluidAccN_;

        //- Secant scale estimates from the FSI iteration pairs of the
        //  current time step
        DynamicList<scalar> secantLambdas_;

        //- Impedance samples (solid impedance, fluid impedance, mean seed
        //  coefficient seen by the sample) of the current time step
        DynamicList<vector> secantSamples_;

        //- Impedance samples of previous time steps and their time indices
        DynamicList<vector> secantSamplesHist_;
        DynamicList<label> secantSamplesHistStep_;

        //- Accumulated normalised per-face least-squares terms of the
        //  current time step
        scalarField secantNumField_;
        scalarField secantDenField_;

        //- Pressure-change norm at the previous FSI iteration
        scalar prevDeltaPNorm_;

        //- Number of successive diverging FSI iterations
        label nDiverging_;

        //- Diagnostics file
        autoPtr<OFstream> diagFilePtr_;


    // Private functions

        //- Return the Robin coefficient rho_s*hs on the fluid patch
        const scalarField& rhoSolidHs() const;

        //- Return the FSI solver
        const fluidSolidInterface& fsi() const;

        //- Return the FSI interface index of this patch
        label interfaceIndex() const;

        //- Map a solid patch face field to this fluid patch
        tmp<scalarField> solidPatchToFluidPatch
        (
            const scalarField& solidField
        ) const;

        //- Return the effective time step of the fluid momentum equation
        scalar effectiveDeltaT() const;

        //- Return the solid wall thickness on the solid patch, found by ray
        //  casting through the solid; wetted-on-both-sides walls are halved
        const scalarField& solidThickness() const;

        //- Calculate the a-priori coefficient rho_s*hs on the fluid patch
        tmp<scalarField> calcSeedCoeff(const word& model) const;

        //- Return the a-priori coefficient, updated once per time step
        const scalarField& seedCoeff() const;

        //- Return the fluid normal acceleration implied by the pressure
        tmp<scalarField> fluidNormalAcceleration() const;

        //- Smooth a face field over the patch face neighbours
        void smoothPatchField(scalarField& fld, const label nSweeps) const;

        //- Bound and apply a new secant scale factor field
        void setSecantScale(const scalarField& newScale);

        //- Update the secant scale from the estimates of the current step
        void applySecantUpdate();

        //- Reset the Robin coefficient state (not the settings)
        void resetRobinState();

        //- Interface fluxes: leakage relative to the mesh motion (sum of
        //  magnitudes and net), mesh motion, and the carry-over mismatch
        //  between the old fluid interface velocity and old mesh velocity
        void leakageFluxes
        (
            scalar& leakFlux,
            scalar& netLeakFlux,
            scalar& wallMotionFlux,
            scalar& carryOverFlux
        ) const;

        //- Return the diagnostics file (master only)
        OFstream& diagFile();

        //- Report the dimensionless Robin number once per time step
        void reportRobinNumber(const scalarField&);


public:

    //- Runtime type information
    TypeName("elasticWallPressure");


    // Constructors

        //- Construct from patch and internal field
        elasticWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        elasticWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given elasticWallPressureFvPatchScalarField
        //  onto a new patch
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new elasticWallPressureFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new elasticWallPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            //- Get previous pressure patch field
            const scalarField& prevPressure() const
            {
                return prevPressure_;
            }

            //- Get previous pressure patch field
            scalarField& prevPressure()
            {
                return prevPressure_;
            }

            //- Get previous pressure patch field
            const vectorField& prevAcceleration() const
            {
                return prevAcceleration_;
            }

            //- Get previous pressure patch field
            vectorField& prevAcceleration()
            {
                return prevAcceleration_;
            }

            //- Return the normal solid acceleration used by the Robin
            //  condition
            tmp<scalarField> prevNormalAcceleration() const
            {
                return patch().nf() & prevAcceleration_;
            }

            //- Set by the FSI coupling: is the fluid interface moved to the
            //  solid interface in every FSI iteration
            void setFluidMeshFollowsSolid(const bool follows)
            {
                fluidMeshFollowsSolid_ = follows;
            }

            //- Is the fluid interface moved to the solid interface in every
            //  FSI iteration
            bool fluidMeshFollowsSolid() const
            {
                return fluidMeshFollowsSolid_;
            }

            //- Set by the fluid model if the interface flux is built such
            //  that the converged flux equals the mesh flux
            void setKinematicConsistency(const bool consistent)
            {
                kinematicConsistency_ = consistent;
            }

            //- Is the interface flux kinematically consistent with the mesh
            //  motion, i.e. does the leakage vanish as the FSI iterations
            //  converge
            bool kinematicConsistency() const
            {
                return kinematicConsistency_;
            }

        // Robin coefficient

            //- Update the Robin coefficient estimate and diagnostics from the
            //  latest solid response. Called by the FSI coupling each outer
            //  iteration after prevPressure and prevAcceleration are set.
            void updateRobinCoefficient();

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate patch flux
            virtual void patchFlux
            (
                GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
                const fvMatrix<scalar>& matrix
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
