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
    Foam::externalPort

Description
    The seam through which a 3D (or any external) domain exchanges flow and
    pressure with the lumped network: one port is one coupled patch.

    The domain supplies the flow Q through its patch and the network answers
    with the pressure to impose there. Q is positive in the natural
    direction of the patch,

        inlet   Q into the domain  (mitral inflow Qmv)
        outlet  Q out of the domain (aortic outflow Qav)

    so the port injects -Q (inlet) or +Q (outlet) into the lumped node it is
    attached to. An optional valve law stands between node and patch, and is
    inverted for the imposed pressure,

        outlet  Ppatch = Pnode + dP(Q)
        inlet   Ppatch = Pnode - dP(Q)

    with dP = orificeValve::dPfromQ, the same bounded bidirectional law the
    left-ventricle coupling uses.

    Where Q comes from is the provider:

        external   set by the coupling code over each step as a polynomial of
                   at most third degree in alpha = (t - t0)/(t1 - t0),
                   Q = c0 + c1 alpha + c2 alpha^2 + c3 alpha^3
        table      read from a waveform table (replay of a recorded run)

        type externalPort;
        direction inlet;                   // or outlet
        patch     outlet;                  // the CFD patch it is coupled to
        provider  external;                // or table
        valve { CQ 400; reverseFraction 0.005; }    // optional
        waveform { type table; file "..."; column 1; }  // provider table

    Ports are normally created by a coupledDomain, not written by hand.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef externalPort_H
#define externalPort_H

#include "circulationComponent.H"
#include "circulationWaveform.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class externalPort Declaration
\*---------------------------------------------------------------------------*/

class externalPort
:
    public circulationComponent
{
public:

    // Public Typedefs

        //- Coefficients c0..c3 of the forcing polynomial in alpha
        typedef FixedList<scalar, 4> forcingCoeffs;


private:

    // Private Data

        //- True for an inlet port, false for an outlet port
        const bool inlet_;

        //- True if the flow is read from a waveform table (provider table)
        const bool tabulated_;

        //- Name of the coupled patch, if given
        const word patch_;

        //- True if a valve law stands between the node and the patch
        bool hasValve_;

        //- Valve coefficient CQ [m3/s/Pa^0.5]
        scalar CQ_;

        //- Regurgitant orifice fraction CQ_r/CQ of the valve
        scalar reverseFraction_;

        //- Flow waveform, for provider table
        autoPtr<circulationWaveform> table_;

        //- Interval [t0, t1] of the external forcing [s]
        scalar t0_, t1_;

        //- External forcing over [t0, t1], coefficients in ascending powers
        //  of alpha
        forcingCoeffs c_;


public:

    //- Runtime type information
    TypeName("externalPort");


    // Constructors

        //- Construct reading direction and provider from the dictionary
        externalPort
        (
            const word& name,
            const dictionary& dict,
            const circulationContext& ctx
        );

        //- Construct with direction and provider given
        externalPort
        (
            const word& name,
            const dictionary& dict,
            const circulationContext& ctx,
            const word& direction,
            const word& provider
        );


    // Member Functions

        //- Role: portSource
        virtual role kind() const { return role::portSource; }

        //- Report the valve coefficient, if there is a valve, for the setup
        //  check
        virtual void parameters(DynamicList<parameter>& p) const
        {
            if (hasValve_)
            {
                p.append
                (
                    {
                        "valve CQ",
                        circulationContext::quantity::valveCoefficient,
                        CQ_
                    }
                );
            }
        }

        //- True for an inlet port
        bool inlet() const { return inlet_; }

        //- Name of the CFD patch this port is coupled to ("patch" entry), or
        //  empty
        const word& patchName() const { return patch_; }

        //- True if the flow is read from a waveform table
        bool tabulated() const { return tabulated_; }

        //- Supply the flow over [t0, t1] (SI, polynomial coefficients)
        void setForcing
        (
            const scalar t0,
            const scalar t1,
            const forcingCoeffs& c
        );

        //- Flow through the patch in its natural direction at time t, from
        //  the table or the forcing polynomial [m3/s]
        virtual scalar portFlow(const scalar t) const;

        //- Time derivative of portFlow, from the forcing polynomial [m3/s2].
        //  At and beyond the ends of [t0, t1] it is the one-sided derivative
        //  there.
        scalar portFlowRate(const scalar t) const;

        //- Flow injected into the node: -portFlow (inlet) or +portFlow
        //  (outlet) [m3/s]
        virtual scalar injectedFlow(const scalar t) const;

        //- Time derivative of injectedFlow [m3/s2]
        virtual scalar injectedFlowRate(const scalar t) const;

        //- Pressure to impose at the patch: Pnode, minus (inlet) or plus
        //  (outlet) the valve pressure drop if there is a valve [Pa]
        virtual scalar patchPressure(const scalar Pnode, const scalar t) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
