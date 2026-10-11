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
    Foam::immersedBodyMotions::fsiDriven

Description
    Deformation of the immersed body given by a fluid-solid interface: the
    positions and velocities of the surface points are set by the interface
    in each coupling iteration (setPointsAndVelocities), and the
    configuration index is incremented each time, so that the body is
    re-positioned in the fluid mesh within the time step. Until they are
    set, the points are at their reference positions and at rest.

    The surface is normally given by the interface too (the body has no
    surface entry); see fluidSolidInterface.

    Usage:
    \verbatim
    motion
    {
        type        fsiDriven;
    }
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    fsiDriven.C

\*---------------------------------------------------------------------------*/

#ifndef fsiDriven_H
#define fsiDriven_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                          Class fsiDriven Declaration
\*---------------------------------------------------------------------------*/

class fsiDriven
:
    public immersedBodyMotion
{
    // Private Data

        //- Current positions of the surface points, empty until set
        pointField points_;

        //- Current velocities of the surface points, empty until set
        vectorField velocities_;

        //- Configuration index, incremented each time the points are set
        label configuration_;


public:

    //- Runtime type information
    TypeName("fsiDriven");


    // Constructors

        //- Construct from dictionary
        explicit fsiDriven(const dictionary& dict);


    //- Destructor
    virtual ~fsiDriven() = default;


    // Member Functions

        //- Set the positions and velocities of the surface points
        void setPointsAndVelocities
        (
            const pointField& points,
            const vectorField& velocities
        );

        //- Positions at time t of the surface points: the positions set by
        //  the interface, else the reference positions
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const;

        //- The motion is not rigid
        virtual bool rigid() const
        {
            return false;
        }

        //- Configuration index
        virtual label configuration() const
        {
            return configuration_;
        }

        //- Not defined for a deforming motion: use pointVelocities
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const;

        //- Velocities of the surface points: those set by the interface,
        //  else zero
        virtual tmp<vectorField> pointVelocities
        (
            const pointField& points0,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBodyMotions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
