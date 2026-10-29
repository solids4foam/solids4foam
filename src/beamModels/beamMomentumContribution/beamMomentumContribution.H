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
    Foam::beamMomentumContribution

Description
    Abstract base class for adding external linear and angular momentum
    contributions to the beam solver.

    This class broadly follows the morisonForce class of Niels Gjoel Jacobsen,
    Deltares.

SourceFiles
    beamMomentumContribution.C

Author
    Philip Cardiff, UCD.
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef beamMomentumContribution_H
#define beamMomentumContribution_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class beamMomentumContribution Declaration
\*---------------------------------------------------------------------------*/

class beamMomentumContribution
// :
// public IOdictionary
{

private:

    // Private data
    const word name_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    beamMomentumContribution(const beamMomentumContribution&);

    //- Disallow default bitwise assignment
    void operator=(const beamMomentumContribution&);

    //- Stop with an error if the beam region holds more than one beam.
    //  Multi-beam support is tracked in issue #46
    void checkSingleBeam(const beamModel& bm) const;

    //- Stop with an error unless the beam has a circular cross-section,
    //  for contributions that use R() as the radius. R() of the other
    //  cross-section models is not a radius (issue #34)
    void checkCircularSection(const beamModel& bm) const;

protected:

    // Protected Member Functions

    //- Return the unit tangents to the beam centreline at the beam
    //  CV cell-centres. They are built from beam 0 alone, so a contribution
    //  that uses them must return true from requiresSingleBeam()
    tmp<vectorField> unitTangents(const beamModel& bm) const;

public:

    //- Runtime type information
    TypeName("beamMomentumContribution");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            beamMomentumContribution,
            dictionary,
            (
                // const Time& runTime
             const word name,
             const dictionary& dict
            ),
            (name, dict)
        );

    // Constructors

        //- Construct from components
        beamMomentumContribution
        (
            const word& name,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<beamMomentumContribution> New
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~beamMomentumContribution();


    // Member Functions

        //- Whether the contribution supports only a single beam per region.
        //  Default: false
        virtual bool requiresSingleBeam() const
        {
            return false;
        }

        //- Whether the contribution uses R() as the radius of a circular
        //  cross-section. Default: false
        virtual bool requiresCircularSection() const
        {
            return false;
        }

        //- Stop with an error if the beam does not suit the contribution,
        //  as declared by requiresSingleBeam() and requiresCircularSection().
        //  Called once by the beam model, after it creates the contributions
        void checkBeam(const beamModel& bm) const;

        //- Hook called once per beam evolve() before the NR loop.
        //  Default no-op.  Contributions that need to perform work that
        //  must remain frozen across the inner NR iterations (e.g. external
        //  fluid sampling) should do it here.
        virtual void preEvolve(const beamModel& bm)
        {}

        //- Return the contribution to the block diagonal of the Jacobian matrix
        //  Each block should be a 6x6 matrix
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        ) = 0;

        //- Return the contribution to the linear momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> linearMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        ) = 0;

        //- Return the contribution to the angular momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> angularMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
