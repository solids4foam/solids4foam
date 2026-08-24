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
    volStrainRateStab

Description
    Face-based stabilisation model targeting temporal oscillations in the
    volumetric strain tr(gradD).

    The stabilisation face vector at face f is:

        faceVector_f = scaleFactor * C(mode) * trDiff_f * n_f

    where

        trDiff_f = tr(gradD_f_future) - tr(gradD_f_old)

    is the face-interpolated difference in volumetric strain between the
    current (future-time extrapolated) and previous time levels, and C(mode)
    is a scaling coefficient that depends on the selected mode:

    Mode                  C(mode)              Vanishes as
    ──────────────────    ─────────────────    ─────────────────────────────────────
    physicalDamping       tau/deltaT           never (permanent physical viscosity)
    firstOrderTemporal    1                    O(deltaT)   as deltaT -> 0
    secondOrderTemporal   deltaT               O(deltaT^2) as deltaT -> 0
    spatioTemporal        h^2/m^2              O(h^2*deltaT) as h,deltaT -> 0
    h2PhysicalDamping     h^2/(m^2*deltaT)     O(h^2) as h->0; O(1) as deltaT->0

    Notes on scaleFactor and dimensional consistency
    ─────────────────────────────────────────────────
    The momentum stabilisation is constructed with dimless (dimensionless)
    dims, so faceVector must be dimensionless.  All modes satisfy this within
    OpenFOAM's type system; however the implicit physical units of scaleFactor
    differ by mode:

      physicalDamping      scaleFactor is dimensionless; tau is in seconds
                           (set to a characteristic damping time scale).
      firstOrderTemporal   scaleFactor is dimensionless.
      secondOrderTemporal  scaleFactor carries an implicit [1/s]; the deltaT
                           factor is extracted as a plain double.
      spatioTemporal       scaleFactor carries an implicit [1/m^2]; h^2 is
                           non-dimensionalised by sqr(dimLength).
      h2PhysicalDamping    scaleFactor is dimensionless; h^2 provides spatial
                           scaling and 1/deltaT provides temporal persistence.

    The scaleFactor (from the base class) and the optional tau are the only
    user-tuning parameters.

    Jacobian approximation
    ──────────────────────
    vectorJacobian provides a Laplacian-based implicit approximation:

        J = scaleFactorJacobian * C(mode) * fvm::laplacian([gamma,] field)

    For modes physicalDamping and secondOrderTemporal, C(mode) depends on
    deltaT; the Jacobian should be rebuilt each time step by passing
    rebuild = true from the calling solid model.

    Dictionary example
    ──────────────────

        type            volStrainRate;
        scaleFactor     0.1;
        mode            firstOrderTemporal;  // or physicalDamping,
                                             //    secondOrderTemporal,
                                             //    spatioTemporal

    For physicalDamping mode, tau is also required:

        type            volStrainRate;
        scaleFactor     0.1;
        mode            physicalDamping;
        tau             1e-4;   // damping time scale in seconds

SourceFiles
    volStrainRateStab.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef volStrainRateStab_H
#define volStrainRateStab_H

#include "stabilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class volStrainRateStab Declaration
\*---------------------------------------------------------------------------*/

class volStrainRateStab
:
    public stabilisationModel
{
    // Private data

        //- Stabilisation mode: physicalDamping | firstOrderTemporal |
        //  secondOrderTemporal | spatioTemporal | h2PhysicalDamping
        const word mode_;

        //- Physical damping time scale [s]; read from dict only when
        //  mode == physicalDamping
        const scalar tau_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        volStrainRateStab(const volStrainRateStab&);

        //- Disallow default bitwise assignment
        void operator=(const volStrainRateStab&);


public:

    //- Runtime type information
    TypeName("volStrainRate");


    // Constructors

        //- Construct from dictionary
        volStrainRateStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~volStrainRateStab();


    // Member Functions

        //- Re-calculate the vector stabilisation.
        //  gradPtr must be non-null; it must point to the gradient of the
        //  displacement field with an accessible old-time level
        //  (gradPtr->oldTime()).
        virtual void updateVector
        (
            const volVectorField& field,
            const volTensorField* gradPtr = nullptr
        ) const;

        //- Provide an implicit (Laplacian-based) approximation of the
        //  stabilisation Jacobian.
        //  An optional face diffusivity gammaPtr may be supplied.
        //  For physicalDamping, secondOrderTemporal and h2PhysicalDamping
        //  modes the coefficient depends on deltaT, so callers should pass
        //  rebuild = true when the time step changes.
        virtual const fvVectorMatrix& vectorJacobian
        (
            const volVectorField& field,
            const surfaceScalarField* gammaPtr = nullptr,
            const bool rebuild = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
