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
    Foam::resistor

Description
    Ohm's law, Q = (Pin - Pout)/R.

        type resistor;
        R    0.5;           // pressure*s/volume

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef resistor_H
#define resistor_H

#include "circulationComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class resistor Declaration
\*---------------------------------------------------------------------------*/

class resistor
:
    public circulationComponent
{
    // Private Data

        //- Resistance R [Pa s/m3]
        scalar R_;


public:

    //- Runtime type information
    TypeName("resistor");


    // Constructors

        //- Construct from dictionary
        resistor
        (
            const word& name,
            const dictionary& dict,
            const circulationContext& ctx
        );


    // Member Functions

        //- Report R for the setup check
        virtual void parameters(DynamicList<parameter>&) const;

        //- Role: branch
        virtual role kind() const { return role::branch; }

        //- The flow is linear in dP
        virtual bool linear() const { return true; }

        //- Flow in -> out, dP/R [m3/s]
        virtual scalar flow(const scalar, const scalar dP) const
        {
            return dP/R_;
        }

        //- d flow / d dP, 1/R [m3/s/Pa]
        virtual scalar dFlowdDP(const scalar, const scalar) const
        {
            return 1/R_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
