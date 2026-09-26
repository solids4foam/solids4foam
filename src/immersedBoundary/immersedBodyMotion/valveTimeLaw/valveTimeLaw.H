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
    Foam::valveTimeLaw

Description
    Periodic gain G(t) in [0, 1], and its time derivative, for the heart
    valve motions, as a function of the phase
    \f$ \phi = (t \bmod T)/T \in [0, 1) \f$ of the cardiac cycle of period T.

    The laws are
    - \c cos2: \f$ G = (1 - \cos 2 \pi \phi)/2 \f$;
    - \c sin2: \f$ G = \sin^2 \pi \phi \f$;
    - \c smoothstep: \f$ G = S(\phi) \f$, with the cubic smoothstep
      \f$ S(u) = 3 u^2 - 2 u^3 \f$;
    - \c twoWindow: \f$ G = E_c (1 - E_o) \f$;
    - \c threeWindow: \f$ G = \min(E_{c1} (1 - E_o) + E_{c2}, 1) \f$;

    where the edge \f$ E_w(\phi) = S((\phi - \phi_0)/(\phi_1 - \phi_0)) \f$
    rises smoothly from 0 to 1 across the window \f$ w = (\phi_0, \phi_1) \f$:
    \c closeWindow (\f$ E_c \f$) or \c close1Window (\f$ E_{c1} \f$),
    \c openWindow (\f$ E_o \f$) and \c close2Window (\f$ E_{c2} \f$), given
    as fractions of the period with \f$ 0 \le \phi_0 \le \phi_1 \le 1 \f$.

    Usage, in the dictionary of the motion:
    \verbatim
    period          0.8;              // [s]
    timeLaw         twoWindow;
    closeWindow     (0.625 0.655);    // twoWindow
    openWindow      (0.94 1);         // twoWindow and threeWindow
    close1Window    (0.1 0.2);        // threeWindow
    close2Window    (0.6 0.7);        // threeWindow
    \endverbatim

    The edges and laws follow the valve time gains of the immersed boundary
    library contributed to cardiacFoam by Sairam Pamulaparthi Venkata; windows
    that wrap around the end of the cycle (\f$ \phi_0 > \phi_1 \f$) are not
    supported.

Author
    Based on the valve time laws of Sairam Pamulaparthi Venkata, UCD.
    Philip Cardiff, UCD.

SourceFiles
    valveTimeLaw.C

\*---------------------------------------------------------------------------*/

#ifndef valveTimeLaw_H
#define valveTimeLaw_H

#include "dictionary.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class valveTimeLaw Declaration
\*---------------------------------------------------------------------------*/

class valveTimeLaw
{
    // Private Data

        //- Period [s]
        scalar period_;

        //- Law
        word law_;

        //- Windows, as fractions of the period
        Pair<scalar> close1_;
        Pair<scalar> open_;
        Pair<scalar> close2_;


public:

    // Static Member Functions

        //- Read the window called name from dict, checking that
        //  0 <= begin <= end <= 1
        static Pair<scalar> readWindow
        (
            const dictionary& dict,
            const word& name
        );

        //- Phase in [0, 1) of time t in a cycle of the given period
        static scalar phase(const scalar t, const scalar period);

        //- Edge rising from 0 to 1 across the window w, at phase phi
        static scalar edge(const scalar phi, const Pair<scalar>& w);

        //- Derivative of the edge with respect to the phase
        static scalar edgeDerivative
        (
            const scalar phi,
            const Pair<scalar>& w
        );


    // Constructors

        //- Construct from the dictionary of the motion
        explicit valveTimeLaw(const dictionary& dict);


    // Member Functions

        //- Period [s]
        scalar period() const
        {
            return period_;
        }

        //- Gain at time t
        scalar value(const scalar t) const;

        //- Time derivative of the gain at time t [1/s]
        scalar derivative(const scalar t) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
