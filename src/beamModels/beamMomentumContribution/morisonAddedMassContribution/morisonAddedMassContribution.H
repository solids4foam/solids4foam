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
    Foam::beamMomentumContribution::morisonAddedMassContribution

Description
    Derived class to evaluate the Morison added-mass force contributions
    which will be added to the beam solver governing eqns

    The beam accelerates through still fluid, so the added-mass force per
    unit length is

        f = -rhoFluid*pi*R^2*(CMn*a_n + CMt*a_t)

    where a_n and a_t are the components of the beam acceleration normal and
    tangential to the beam centreline, and CMn and CMt are the normal and
    tangential added-mass coefficients (e.g. CMn = 1, CMt = 0 for a mooring
    chain). The acceleration is evaluated as in the beam inertia term, so the
    force is added to the linear momentum source, and its linearisation to
    the diagonal of the Jacobian, in the same way as the beam inertia.

    Notes:
    - CMn and CMt are the added-mass coefficients C_a (as Can and Cat in
      MoorDyn), not the Morison inertia coefficient C_M = 1 + C_a.
    - Only the beam's own acceleration is used. The fluid-acceleration
      (Froude-Krylov) terms rhoFluid*V*(1 + C_a)*a_fluid are not included,
      so in waves, or in a coupled free-surface simulation, the contribution
      is incomplete.
    - rhoFluid applies along the whole beam, including any part above a free
      surface.
    - Only valid for a circular cross-section (crossSectionModel circle):
      pi*R^2 is the displaced area and the added mass is the same in all
      directions normal to the beam. Other cross-sections stop with an error.
    - A beam region must hold a single beam.

    Example of an entry in constant/beamMomentumContributionProperties:
    \verbatim
    morisonAddedMass
    {
        beamMomentumContributionType morisonAddedMass;
        morisonAddedMassCoeffs
        {
            CMn       1.0;      // normal added-mass coefficient
            CMt       0.0;      // tangential added-mass coefficient
            rhoFluid  1000;     // fluid density [kg/m^3]
        }
    }
    \endverbatim

SourceFiles
    morisonAddedMassContribution.C

Author
    Amirhossein Taran, UCD.


\*---------------------------------------------------------------------------*/

#ifndef morisonAddedMassContribution_H
#define morisonAddedMassContribution_H

#include "beamMomentumContribution.H"
#include "HermiteSpline.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class morisonAddedMassContribution Declaration
\*---------------------------------------------------------------------------*/

class morisonAddedMassContribution
:
    public beamMomentumContribution
{

private:

    // Private Data

        //- beam momentum contributions sub dictionary
        dictionary beamMomentumContribDict_;

        //- Normal added-mass coefficient for Morison Force
        scalar CMn_;

        //- Tangential added-mass coefficient for Morison Force
        scalar CMt_;

        //- Fluid density
        scalar rhoFluid_;


    // Private Member Functions

        //- Return the d2dt2 scheme of W, as selected by the beam model.
        //  Stops with an error unless it is Euler or Newmark
        static word d2dt2SchemeName(const fvMesh& mesh);

        //- Return the displaced fluid mass of each beam CV,
        //  rhoFluid*pi*R^2*L. Stops with an error unless the beam has a
        //  circular cross-section
        tmp<scalarField> fluidMass(const beamModel& bm) const;


public:

    //- Runtime type information
    TypeName("morisonAddedMass");

    // Constructors

        //- Construct from components
        morisonAddedMassContribution
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        ~morisonAddedMassContribution()
        {}


    // Member Functions

        //- The tangents and the radius come from beam 0 alone
        virtual bool requiresSingleBeam() const
        {
            return true;
        }

        //- R() is used as the radius
        virtual bool requiresCircularSection() const
        {
            return true;
        }

        //- Return the contribution to the block diagonal of the Jacobian matrix
        //  Each block should be a 6x6 matrix
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        // bm: Reference to the beamModel object
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );

        //- Return the contribution to the linear momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        // bm: Reference to the beamModel object
        virtual tmp<vectorField> linearMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );

        //- Return the contribution to the angular momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        // bm: Reference to the beamModel object
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
