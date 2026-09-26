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
    Foam::immersedBodyMotion

Description
    Run-time selectable prescribed motion of an immersed body surface.

    A motion gives the position of the surface points at a given time, from
    their reference (as read) positions, and the velocity of the body at
    given points in the current configuration, e.g. the cell centres covered
    by the body.

Author
    Philip Cardiff, UCD.

SourceFiles
    immersedBodyMotion.C
    immersedBodyMotionNew.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBodyMotion_H
#define immersedBodyMotion_H

#include "dictionary.H"
#include "pointField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class immersedBodyMotion Declaration
\*---------------------------------------------------------------------------*/

class immersedBodyMotion
{
protected:

    // Protected Data

        //- Motion coefficients
        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("immersedBodyMotion");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            immersedBodyMotion,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from dictionary
        explicit immersedBodyMotion(const dictionary& dict);

        //- No copy construct
        immersedBodyMotion(const immersedBodyMotion&) = delete;

        //- No copy assignment
        void operator=(const immersedBodyMotion&) = delete;


    // Selectors

        //- Select the motion given by the "type" entry of dict
        static autoPtr<immersedBodyMotion> New(const dictionary& dict);


    //- Destructor
    virtual ~immersedBodyMotion() = default;


    // Member Functions

        //- Does the body move? A body that does not move is placed once
        virtual bool moving() const
        {
            return true;
        }

        //- Is the motion rigid? The velocity of a rigid motion is defined
        //  at any point; that of a deforming motion only at the surface
        //  points, from pointVelocities
        virtual bool rigid() const
        {
            return true;
        }

        //- Positions at time t of the surface points whose reference
        //  positions are points0
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const = 0;

        //- Velocity of the body at time t at the points x, given in the
        //  current configuration (rigid motions)
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const = 0;

        //- Velocity at time t of the surface points whose reference
        //  positions are points0. The default is the velocity at their
        //  positions at time t
        virtual tmp<vectorField> pointVelocities
        (
            const pointField& points0,
            const scalar t
        ) const
        {
            return velocity(points(points0, t), t);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
