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
    Foam::beamMomentumContribution::morisonDragContribution

Description
    Derived class to evaluate the (explicit) Morison drag force
    contributions which will be added to the beam solver governing eqns

    Only valid for a circular cross-section (crossSectionModel circle),
    as R is used as the radius. A beam region must hold a single beam.
    Other cases stop with an error.

SourceFiles
    morisonDragContribution.C

Author
    Seevani Bali, UCD.
    Philip Cardiff, UCD.
    Amirhossein Taran, UCD.


\*---------------------------------------------------------------------------*/

#ifndef morisonDragContribution_H
#define morisonDragContribution_H

#include "beamMomentumContribution.H"
#include "HermiteSpline.H"
// #include "beamModel.H"
#include "coupledTotalLagNewtonRaphsonBeam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// namespace beamMomentumContributions
// {

/*---------------------------------------------------------------------------*\
                           Class morisonDragContribution Declaration
\*---------------------------------------------------------------------------*/

class morisonDragContribution
:
    public beamMomentumContribution
{

private:

    //- beam momentum contributions sub dictionary
    dictionary beamMomentumContribDict_;

    //- Normal drag coefficient for Morison Force
    scalar Cdn_;

    //- Tangential drag coefficient for Morison Force
    scalar Cdt_;

    //- Fluid density
    scalar rhoFluid_;

public:

    //- Runtime type information
    TypeName("morisonDrag");

    // Constructors

        //- Contruct from components
        morisonDragContribution
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        ~morisonDragContribution()
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

// } // End namespace beamMomentumContributions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
