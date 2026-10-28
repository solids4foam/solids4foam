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

Description
    The flow, the injected flow and the imposed patch pressure of an
    external coupling port, including the valve law between the lumped node
    and the patch, and the forcing set by the coupling code.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "externalPort.H"
#include "orificeValve.H"

namespace Foam
{
    defineTypeNameAndDebug(externalPort, 0);
    addToRunTimeSelectionTable
    (
        circulationComponent,
        externalPort,
        dictionary
    );
}

namespace
{
    bool readInlet(const Foam::word& direction, const Foam::dictionary& d)
    {
        if (direction == "inlet")
        {
            return true;
        }

        if (direction == "outlet")
        {
            return false;
        }

        FatalIOErrorInFunction(d)
            << "direction must be inlet or outlet, not '" << direction << "'"
            << Foam::exit(Foam::FatalIOError);
        return false;
    }

    bool readTabulated(const Foam::word& provider, const Foam::dictionary& d)
    {
        if (provider == "table")
        {
            return true;
        }

        if (provider == "external")
        {
            return false;
        }

        FatalIOErrorInFunction(d)
            << "provider must be external or table, not '" << provider << "'"
            << Foam::exit(Foam::FatalIOError);
        return false;
    }
}

Foam::externalPort::externalPort
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    externalPort
    (
        name,
        dict,
        ctx,
        dict.get<word>("direction"),
        dict.getOrDefault<word>("provider", "external")
    )
{}


Foam::externalPort::externalPort
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx,
    const word& direction,
    const word& provider
)
:
    circulationComponent(name, dict, ctx),
    inlet_(readInlet(direction, dict)),
    tabulated_(readTabulated(provider, dict)),
    patch_(dict.getOrDefault<word>("patch", word::null)),
    hasValve_(dict.isDict("valve")),
    CQ_(0),
    reverseFraction_(0.005),
    t0_(0),
    t1_(1),
    c_(0.0)
{
    if (hasValve_)
    {
        const dictionary& v = dict.subDict("valve");
        CQ_ = ctx.read(v, "CQ", circulationContext::quantity::valveCoefficient);
        reverseFraction_ = v.getOrDefault<scalar>("reverseFraction", 0.005);

        if (CQ_ <= 0)
        {
            fail("valve CQ must be positive");
        }
    }

    if (tabulated_)
    {
        table_.reset
        (
            new circulationWaveform
            (
                dict,
                ctx,
                circulationContext::quantity::flow
            )
        );
    }
}


void Foam::externalPort::setForcing
(
    const scalar t0,
    const scalar t1,
    const forcingCoeffs& c
)
{
    if (tabulated_)
    {
        fail("is table driven; its flow cannot be set from outside");
    }

    t0_ = t0;
    t1_ = t1;
    c_ = c;
}


Foam::scalar Foam::externalPort::portFlow(const scalar t) const
{
    if (tabulated_)
    {
        return table_().value(t);
    }

    const scalar dt = t1_ - t0_;
    const scalar a =
        dt > VSMALL
      ? Foam::min(Foam::max((t - t0_)/dt, scalar(0)), scalar(1))
      : scalar(0);

    return c_[0] + a*(c_[1] + a*(c_[2] + a*c_[3]));
}


Foam::scalar Foam::externalPort::injectedFlow(const scalar t) const
{
    const scalar Q = portFlow(t);
    return inlet_ ? -Q : Q;
}


Foam::scalar Foam::externalPort::patchPressure
(
    const scalar Pnode,
    const scalar t
) const
{
    if (!hasValve_)
    {
        return Pnode;
    }

    const scalar dP = orificeValve::dPfromQ(CQ_, portFlow(t), reverseFraction_);

    return inlet_ ? Pnode - dP : Pnode + dP;
}

// ************************************************************************* //
