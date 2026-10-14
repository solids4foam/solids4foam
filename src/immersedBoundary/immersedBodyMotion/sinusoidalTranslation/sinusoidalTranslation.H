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
    Foam::immersedBodyMotions::sinusoidalTranslation

Description
    Sinusoidal translation of the immersed body:

    \f[
        \vec{x}(t) = \vec{x}_0 + A \sin(\omega t + \phi) \vec{d} + \vec{s}
    \f]

    with velocity \f$ A \omega \cos(\omega t + \phi) \vec{d} \f$, where
    \f$ \omega = 2 \pi / T \f$, \f$ T \f$ is the period, \f$ \phi \f$ the
    phase, \f$ \vec{d} \f$ the (normalised) direction and \f$ \vec{s} \f$ a
    constant offset.

    Usage:
    \verbatim
    motion
    {
        type        sinusoidalTranslation;
        amplitude   0.25;      // [m]
        period      4;         // [s]
        phase       0;         // [rad], optional, default 0
        direction   (1 0 0);   // optional, default (1 0 0)
        offset      (0 0 0);   // [m], optional, default (0 0 0)
    }
    \endverbatim

Author
    Based on the sinusoidal translation of the immersed boundary library
    contributed to cardiacFoam by Sairam Pamulaparthi Venkata, UCD.
    Philip Cardiff, UCD.

SourceFiles
    sinusoidalTranslation.C

\*---------------------------------------------------------------------------*/

#ifndef sinusoidalTranslation_H
#define sinusoidalTranslation_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                   Class sinusoidalTranslation Declaration
\*---------------------------------------------------------------------------*/

class sinusoidalTranslation
:
    public immersedBodyMotion
{
    // Private Data

        //- Amplitude [m]
        const scalar amplitude_;

        //- Period [s]
        const scalar period_;

        //- Phase [rad]
        const scalar phase_;

        //- Direction (normalised)
        const vector direction_;

        //- Constant offset [m]
        const vector offset_;


    // Private Member Functions

        //- Angular frequency
        scalar omega() const;


public:

    //- Runtime type information
    TypeName("sinusoidalTranslation");


    // Constructors

        //- Construct from dictionary
        explicit sinusoidalTranslation(const dictionary& dict);


    //- Destructor
    virtual ~sinusoidalTranslation() = default;


    // Member Functions

        //- Positions at time t of the surface points
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const;

        //- Velocity of the body at time t at the points x
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBodyMotions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
