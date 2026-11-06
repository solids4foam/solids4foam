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
    Foam::compliance

Description
    A compliant vessel: a node that stores volume, P = (V - V0)/C + P0. The
    volume is a state, dV/dt = flow in - flow out.

        type compliance;
        C    1.6;           // volume/pressure
        V0   0;             // unstressed volume, optional (default 0)
        P0   0;             // pressure offset, optional (default 0)
        initialPressure 100;  // or initialVolume; default is V = V0

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef compliance_H
#define compliance_H

#include "circulationComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class compliance Declaration
\*---------------------------------------------------------------------------*/

class compliance
:
    public circulationComponent
{
    // Private Data

        //- Compliance C [m3/Pa]
        scalar C_;

        //- Unstressed volume V0 [m3]
        scalar V0_;

        //- Pressure offset P0 [Pa]
        scalar P0_;

        //- Initial volume [m3]
        scalar Vinit_;


public:

    //- Runtime type information
    TypeName("compliance");


    // Constructors

        //- Construct from dictionary
        compliance
        (
            const word& name,
            const dictionary& dict,
            const circulationContext& ctx
        );


    // Member Functions

        //- Report C, V0 and P0 for the setup check
        virtual void parameters(DynamicList<parameter>&) const;

        //- Role: storage
        virtual role kind() const { return role::storage; }

        //- Pressure for volume V, (V - V0)/C + P0 [Pa]
        virtual scalar pressure(const scalar, const scalar V) const
        {
            return (V - V0_)/C_ + P0_;
        }

        //- Initial volume [m3]
        virtual scalar initialVolume() const { return Vinit_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
