/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::orificeValve

Description
    An orifice valve (the thesis HeartValve): forward flow Q = CQ sqrt(dP)
    for dP = Pin - Pout > 0, no flow otherwise.

        type orificeValve;
        CQ   350;                // volume/(s sqrt(pressure)), or give
        Rhat 5.3e9;              // dP = Rhat Q^2 instead (CQ = 1/sqrt(Rhat))
        R    0;                  // optional resistance in series: dP = Rhat Q^2 + R Q
        epsilon 0;               // optional smoothing pressure; 0 = hard diode
        reverseFraction 0;       // optional regurgitant orifice CQ_r/CQ

    epsilon > 0 replaces the hard switch by a smooth one (needed by implicit
    ODE solvers). reverseFraction > 0 gives a closed valve a small leak,
    Q = -reverseFraction CQ sqrt(-dP).

    The static dPfromQ() is the closed-form inverse used wherever the flow is
    known and the pressure drop is wanted, such as a coupled 3D domain:

        dP =  (Q/CQ)^2                       Q >= 0
        dP = -(|Q|/(fr CQ))^2                Q <  0,   fr in [1e-3, 1]

\*---------------------------------------------------------------------------*/

#ifndef orificeValve_H
#define orificeValve_H

#include "circulationComponent.H"

namespace Foam
{

class orificeValve
:
    public circulationComponent
{
    scalar CQ_;
    scalar Rs_;                 //!< series resistance [Pa s/m3]
    scalar eps_;
    scalar reverseFraction_;

public:

    TypeName("orificeValve");

    orificeValve(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::branch; }
    virtual scalar flow(const scalar t, const scalar dP) const;
    virtual scalar dFlowdDP(const scalar t, const scalar dP) const;

    //- Pressure drop that carries flow Q through the valve [Pa]
    static scalar dPfromQ(const scalar CQ, const scalar Q, const scalar fr);
};

} // End namespace Foam

#endif
