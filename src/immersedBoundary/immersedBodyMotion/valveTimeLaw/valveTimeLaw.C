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

\*---------------------------------------------------------------------------*/

#include "valveTimeLaw.H"
#include "mathematicalConstants.H"
#include "MinMax.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

Foam::Pair<Foam::scalar> Foam::valveTimeLaw::readWindow
(
    const dictionary& dict,
    const word& name
)
{
    const Pair<scalar> w(dict.get<Pair<scalar>>(name));

    if (w.first() < 0 || w.first() > w.second() || w.second() > 1)
    {
        FatalIOErrorInFunction(dict)
            << "The window " << name << " " << w << " must satisfy "
            << "0 <= begin <= end <= 1, as fractions of the period"
            << exit(FatalIOError);
    }

    return w;
}


Foam::scalar Foam::valveTimeLaw::phase(const scalar t, const scalar period)
{
    const scalar phi = t/period - std::floor(t/period);

    return min(max(phi, 0), 1);
}


Foam::scalar Foam::valveTimeLaw::edge
(
    const scalar phi,
    const Pair<scalar>& w
)
{
    if (phi <= w.first())
    {
        return 0;
    }
    else if (phi >= w.second())
    {
        return 1;
    }

    const scalar u = (phi - w.first())/(w.second() - w.first());

    return u*u*(3 - 2*u);
}


Foam::scalar Foam::valveTimeLaw::edgeDerivative
(
    const scalar phi,
    const Pair<scalar>& w
)
{
    if (phi <= w.first() || phi >= w.second())
    {
        return 0;
    }

    const scalar len = w.second() - w.first();
    const scalar u = (phi - w.first())/len;

    return 6*u*(1 - u)/len;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::valveTimeLaw::valveTimeLaw(const dictionary& dict)
:
    period_(dict.getCheck<scalar>("period", scalarMinMax::ge(SMALL))),
    law_(dict.get<word>("timeLaw")),
    close1_(0, 0),
    open_(0, 0),
    close2_(0, 0)
{
    if (law_ == "twoWindow")
    {
        close1_ = readWindow(dict, "closeWindow");
        open_ = readWindow(dict, "openWindow");
    }
    else if (law_ == "threeWindow")
    {
        close1_ = readWindow(dict, "close1Window");
        open_ = readWindow(dict, "openWindow");
        close2_ = readWindow(dict, "close2Window");
    }
    else if (law_ != "cos2" && law_ != "sin2" && law_ != "smoothstep")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown timeLaw " << law_ << nl
            << "Valid laws are: cos2 sin2 smoothstep twoWindow threeWindow"
            << exit(FatalIOError);
    }

    // The gain jumps at the end of each cycle if it does not return to its
    // value at the start
    const scalar G0 = value(0);
    const scalar G1 = value((1 - SMALL)*period_);
    if (mag(G1 - G0) > 1e-6)
    {
        WarningInFunction
            << "The valve time law " << law_ << " goes from " << G0
            << " at the start of the cycle to " << G1 << " at its end: the "
            << "valve jumps between these states at the end of each cycle"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::valveTimeLaw::value(const scalar t) const
{
    const scalar phi = phase(t, period_);

    if (law_ == "twoWindow")
    {
        return edge(phi, close1_)*(1 - edge(phi, open_));
    }
    else if (law_ == "threeWindow")
    {
        return
            min
            (
                edge(phi, close1_)*(1 - edge(phi, open_))
              + edge(phi, close2_),
                1
            );
    }
    else if (law_ == "smoothstep")
    {
        return phi*phi*(3 - 2*phi);
    }

    // cos2 and sin2 are the same law: (1 - cos(2 pi phi))/2 = sin^2(pi phi)
    return sqr(Foam::sin(constant::mathematical::pi*phi));
}


Foam::scalar Foam::valveTimeLaw::derivative(const scalar t) const
{
    const scalar phi = phase(t, period_);

    scalar dGdPhi = 0;

    if (law_ == "twoWindow")
    {
        dGdPhi =
            edgeDerivative(phi, close1_)*(1 - edge(phi, open_))
          - edge(phi, close1_)*edgeDerivative(phi, open_);
    }
    else if (law_ == "threeWindow")
    {
        const scalar G =
            edge(phi, close1_)*(1 - edge(phi, open_)) + edge(phi, close2_);

        if (G < 1)
        {
            dGdPhi =
                edgeDerivative(phi, close1_)*(1 - edge(phi, open_))
              - edge(phi, close1_)*edgeDerivative(phi, open_)
              + edgeDerivative(phi, close2_);
        }
    }
    else if (law_ == "smoothstep")
    {
        dGdPhi = 6*phi*(1 - phi);
    }
    else
    {
        dGdPhi =
            constant::mathematical::pi
           *Foam::sin(constant::mathematical::twoPi*phi);
    }

    return dGdPhi/period_;
}


// ************************************************************************* //
