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
    Foam::elastanceChamber

Description
    A heart chamber with time-varying elastance (Suga/Sagawa),

        P = (V - V0) E(t) + P0,    E = Emin + (Emax - Emin) f(t)

    with the half-cosine activation of the thesis HeartChamber,

        ti = mod(t + (1 - Eshift) T, T)
        f  = (1 - cos(pi ti/Tes))/2                      ti <= Tes
           = (1 + cos(pi (ti - Tes)/(Ted - Tes)))/2      Tes < ti <= Ted
           = 0                                           ti > Ted

        type elastanceChamber;
        V0 5; P0 1; Emin 0.1; Emax 2.5;     // volume, pressure, pressure/volume
        Tes 0.3; Ted 0.45;                  // seconds
        Eshift 0;                           // seconds, optional
        T 1;                                // optional, default Tcycle
        initialVolume 500;

    The atria are the same law with Tes = Tpww/2, Ted = Tpww and
    Eshift = Tpwb.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef elastanceChamber_H
#define elastanceChamber_H

#include "circulationComponent.H"

namespace Foam
{

class elastanceChamber
:
    public circulationComponent
{
    scalar V0_, P0_, Emin_, Emax_;
    scalar Tes_, Ted_, Eshift_, T_;
    scalar Vinit_;

public:

    TypeName("elastanceChamber");

    elastanceChamber
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::storage; }

    //- Elastance at time t [Pa/m3]
    scalar elastance(const scalar t) const;

    virtual scalar pressure(const scalar t, const scalar V) const
    {
        return (V - V0_)*elastance(t) + P0_;
    }

    virtual scalar initialVolume() const { return Vinit_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
