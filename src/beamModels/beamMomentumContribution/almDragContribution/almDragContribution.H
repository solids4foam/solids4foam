/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::almDragContribution

Description
    ALM-style drag contribution for the beam solver.

    Samples the fluid velocity at upstream points (one per beam cell),
    computes a quadratic drag line force per unit length using normal and
    tangential drag coefficients, and assembles it into the beam linear
    momentum source.

    Full form (per beam cell, force per unit length):

        Urel = Uf - Ub                          relative velocity
        Ut   = (Urel . that) that               tangential component
        Un   = Urel - Ut                        normal component

        fn   = -1/2 rho Cdn D |Urel| Un         normal drag line force
        ft   = -1/2 rho Cdt D |Urel| Ut         tangential drag line force

        f    = fn + ft                          total drag line force

    where:
        Uf   = sampled fluid velocity at the upstream point
        Ub   = beam cell velocity
        that = unit beam tangent (dR/ds normalised)
        rho  = fluid density (bm.rhoFluid())
        Cdn  = normal drag coefficient (Cdn_)
        Cdt  = tangential drag coefficient (Cdt_)
        D    = drag reference length (dragReferenceLength_, default 2*bm.R())

    The defaults that use bm.R() (D, and the sampling reference length)
    need a circular cross-section; for other sections, set both
    almSamplingReferenceLength and almDragReferenceLength.

    The stored line force is under-relaxed against the previous estimate:

        almForce = w f + (1 - w) almForce_prev

    with w = forceRelaxation_.  Multiplying almForce by the cell length L
    (in linearMomentumSource) gives the nodal force assembled into the beam
    linear momentum source.

    The sampling is performed once per beam evolve() call (in preEvolve)
    and the per-cell line force is then frozen for the duration of the
    inner NR loop.  Sampling at each NR iteration would cause the upstream
    sample point to wander across fluid cell boundaries and prevent
    convergence (a perfect two-cycle).  Under-relaxation between successive
    outer-corrector forces is governed by forceRelaxation.

    The resulting per-cell line force is written into the volVectorField
    "almForce" registered on the beam mesh, so that downstream consumers
    (e.g. moorFV's beamActuatorLine fvOption) keep working unchanged.

SourceFiles
    almDragContribution.C

Author
    Amirhossein Taran, UCD.

\*---------------------------------------------------------------------------*/

#ifndef almDragContribution_H
#define almDragContribution_H

#include "beamMomentumContribution.H"
#include "meshSearch.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class almDragContribution Declaration
\*---------------------------------------------------------------------------*/

class almDragContribution
:
    public beamMomentumContribution
{

private:

    //- Sub-dictionary of coefficients (name + "Coeffs")
    dictionary coeffs_;

    //- Normal drag coefficient
    scalar Cdn_;

    //- Tangential drag coefficient
    scalar Cdt_;

    //- Under-relaxation factor applied between successive outer-corrector
    //  line force estimates
    scalar forceRelaxation_;

    //- Sampling radius multiplier (passed to getFluidVelocity)
    scalar samplingRadius_;

    //- Reference length used to position the upstream sampling point.
    //  Negative => default to bm.R().
    scalar samplingReferenceLength_;

    //- Reference length used in the drag formula (typically diameter).
    //  Negative => default to 2*bm.R().
    scalar dragReferenceLength_;

    //- Sampling plane normal (must be non-zero)
    vector samplingPlaneNormal_;

    //- Global upstream direction (direction opposite to mean flow)
    vector upstreamDir_;

    //- Whether ground contact constrains the sampling (passed through to
    //  getFluidVelocity)
    bool groundContactActive_;

    //- Ground plane used by the sampling routine: groundPoint and
    //  groundNormal, or groundZ for the plane z = groundZ (default 0)
    plane groundPlane_;

    //- Lazily constructed fluid mesh search engine
    mutable autoPtr<meshSearch> searchEnginePtr_;

    //- Has the decomposition check been run? It is run once, on the first
    //  preEvolve, where the beam mesh is available
    bool decompositionChecked_;

    //- Disallow default bitwise copy construct
    almDragContribution(const almDragContribution&);

    //- Disallow default bitwise assignment
    void operator=(const almDragContribution&);


public:

    //- Runtime type information
    TypeName("almDrag");

    // Constructors

        //- Construct from components
        almDragContribution
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~almDragContribution()
        {}


    // Member Functions

        //- The tangents come from beam 0 alone
        virtual bool requiresSingleBeam() const
        {
            return true;
        }

        //- R() is used as the radius when a reference length is not given
        virtual bool requiresCircularSection() const
        {
            return samplingReferenceLength_ <= 0 || dragReferenceLength_ <= 0;
        }

        //- Sample the fluid and update the cached almForce field
        virtual void preEvolve(const beamModel& bm);

        //- Diagonal Jacobian contribution (currently zero - drag is explicit)
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );

        //- Linear momentum source: returns almForce * L per cell
        virtual tmp<vectorField> linearMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );

        //- Angular momentum source: zero
        virtual tmp<vectorField> angularMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
