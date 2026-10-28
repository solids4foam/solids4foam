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
    Foam::mixedChamber

Description
    A heart chamber whose pressure is an active contraction term switched by
    an activation function plus a passive exponential filling term,

        dV = V - V0
        P  = A(t) (E1 dV + E2 dV^2) + P0 (exp(K dV) - 1)

    (the form of Mynard and of ModularCirc's mixed chambers, and of the
    closed-loop pipe model). Two activations:

        activation { type sineSquared; offset 0.544; duration 0.24; }
            A = sin^2(pi phase/duration) for phase = mod(t - offset, T)
            in [0, duration], zero otherwise
        activation { type halfCosine; Tes 0.3; Ted 0.45; Eshift 0; }
            the half-cosine of elastanceChamber, 0..1..0

    T is the cardiac period (Tcycle, or "T" in the activation dictionary).

        type mixedChamber;
        V0   4e-6;                          // unstressed volume
        active  { E1 2.467e9; E2 -5.6e12; } // pressure/volume,
                                            // pressure/volume^2
        passive { P0 120; K 6.2e4; }        // pressure, 1/volume
        activation { type sineSquared; offset 0; duration 0.4; }
        initialVolume 7.9e-5;

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mixedChamber_H
#define mixedChamber_H

#include "circulationComponent.H"

namespace Foam
{

class mixedChamber
:
    public circulationComponent
{
    scalar V0_, E1_, E2_, P0_, K_;
    scalar Vinit_;

    bool sineSquared_;
    scalar offset_, duration_, T_;
    scalar Tes_, Ted_, Eshift_;

public:

    TypeName("mixedChamber");

    mixedChamber
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::storage; }

    //- Activation 0..1 at time t
    scalar activation(const scalar t) const;

    virtual scalar pressure(const scalar t, const scalar V) const;

    virtual scalar initialVolume() const { return Vinit_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
