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
    Assembles the circulation network from its dictionary and solves it:
    classifies the nets, solves the algebraic nodal pressures, provides the
    derivatives and Jacobian of the ODE system, advances it in time, and
    handles the port forcing, the restart state and the set-up check.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationNetwork.H"

#include <numeric>
#include <sstream>

// * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- Solve A x = b in place (b becomes x), partial pivoting. Returns false if
//  the matrix is numerically singular.
bool solveDense(scalarSquareMatrix& A, scalarField& b)
{
    const label n = b.size();

    scalar scale = 0;
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            scale = Foam::max(scale, Foam::mag(A(i, j)));
        }
    }

    for (label k = 0; k < n; ++k)
    {
        label p = k;
        for (label i = k + 1; i < n; ++i)
        {
            if (Foam::mag(A(i, k)) > Foam::mag(A(p, k)))
            {
                p = i;
            }
        }

        if (Foam::mag(A(p, k)) <= 1e-14*scale || scale == 0)
        {
            return false;
        }

        if (p != k)
        {
            for (label j = 0; j < n; ++j)
            {
                Foam::Swap(A(k, j), A(p, j));
            }
            Foam::Swap(b[k], b[p]);
        }

        for (label i = k + 1; i < n; ++i)
        {
            const scalar f = A(i, k)/A(k, k);
            if (f == 0)
            {
                continue;
            }
            for (label j = k; j < n; ++j)
            {
                A(i, j) -= f*A(k, j);
            }
            b[i] -= f*b[k];
        }
    }

    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label j = i + 1; j < n; ++j)
        {
            s -= A(i, j)*b[j];
        }
        b[i] = s/A(i, i);
    }

    return true;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationNetwork::circulationNetwork
(
    const dictionary& dict,
    const scalar t0
)
:
    ctx_(dict),
    nNets_(0),
    nAlg_(0),
    algLinear_(true),
    nStates_(0),
    t_(t0),
    dxEst_(1e-4)
{
    build(dict);
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::circulationNetwork::build(const dictionary& dict)
{
    // ---- 1. Components -----------------------------------------------------

    circulationBuilder builder(ctx_);

    const dictionary& compsDict = dict.subDict("components");

    for (const entry& e : compsDict)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(compsDict)
                << "Component '" << e.keyword() << "' must be a dictionary"
                << exit(FatalIOError);
        }

        if (e.keyword().find('.') != std::string::npos)
        {
            FatalIOErrorInFunction(compsDict)
                << "Component name '" << e.keyword() << "' must not contain '.'"
                << exit(FatalIOError);
        }

        builder.add
        (
            circulationComponent::New(e.keyword(), e.dict(), ctx_)
        );
    }

    PtrList<circulationComponent>& prims = builder.primitives();

    const label nComps = prims.size();

    if (nComps == 0)
    {
        FatalIOErrorInFunction(dict) << "No components defined"
            << exit(FatalIOError);
    }

    // ---- 2. Terminals and nets ---------------------------------------------

    HashTable<label> terminal;      // "<comp>.in|out" -> terminal id
    forAll(prims, i)
    {
        terminal.insert(prims[i].name() + ".in", 2*i);
        terminal.insert(prims[i].name() + ".out", 2*i + 1);
    }

    labelList parent(2*nComps);
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](label a)
    {
        while (parent[a] != a)
        {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };

    auto unite = [&](const label a, const label b)
    {
        parent[find(a)] = find(b);
    };

    auto lookup = [&](const word& name) -> label
    {
        const word resolved = builder.resolve(name);
        const auto iter = terminal.cfind(resolved);

        if (!iter.good())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown terminal '" << name << "'"
                << (resolved != name ? " (-> " + resolved + ")" : string(""))
                << ". Terminals are <component>.in and <component>.out;"
                << " known components: " << nl
                << wordList(compsDict.toc()) << nl
                << exit(FatalIOError);
        }
        return iter.val();
    };

    // Connections inside composites
    for (const Pair<word>& link : builder.links())
    {
        unite(lookup(link.first()), lookup(link.second()));
    }

    // Connections of the case: chains of components or terminals
    if (dict.found("connections"))
    {
        List<wordList> chains(dict.lookup("connections"));

        for (const wordList& chain : chains)
        {
            if (chain.size() < 2)
            {
                FatalIOErrorInFunction(dict)
                    << "A connection needs at least two names, got "
                    << chain << exit(FatalIOError);
            }

            for (label i = 0; i + 1 < chain.size(); ++i)
            {
                const word& a = chain[i];
                const word& b = chain[i + 1];

                const bool aExplicit = a.find('.') != std::string::npos;
                const bool bExplicit = b.find('.') != std::string::npos;

                unite
                (
                    lookup(aExplicit ? a : word(a + ".out")),
                    lookup(bExplicit ? b : word(b + ".in"))
                );
            }
        }
    }

    // Single-node elements have one node
    forAll(prims, i)
    {
        if (prims[i].singleNode())
        {
            unite(2*i, 2*i + 1);
        }
    }

    labelList netOfRoot(2*nComps, -1);
    nNets_ = 0;
    labelList netOfTerminal(2*nComps);
    forAll(netOfTerminal, t)
    {
        const label r = find(t);
        if (netOfRoot[r] < 0)
        {
            netOfRoot[r] = nNets_++;
        }
        netOfTerminal[t] = netOfRoot[r];
    }

    netIn_.setSize(nComps);
    netOut_.setSize(nComps);
    roles_.setSize(nComps);
    forAll(prims, i)
    {
        netIn_[i] = netOfTerminal[2*i];
        netOut_[i] = netOfTerminal[2*i + 1];
        roles_[i] = prims[i].kind();
    }

    comps_.transfer(prims);

    // ---- 3. Classify nets -------------------------------------------------

    netStorage_.setSize(nNets_, -1);
    netFixed_.setSize(nNets_, -1);
    algIndex_.setSize(nNets_, -1);

    labelList attachments(nNets_, 0);

    forAll(comps_, k)
    {
        const label n = netIn_[k];

        attachments[n]++;
        if (netOut_[k] != n)
        {
            attachments[netOut_[k]]++;
        }

        if (roles_[k] == role::storage)
        {
            if (netStorage_[n] >= 0 || netFixed_[n] >= 0)
            {
                FatalIOErrorInFunction(dict)
                    << "'" << comps_[k].name() << "' shares a node with another"
                    << " compliance, chamber or pressure source. A node can"
                    << " hold only one. Members: " << netMembers(n)
                    << exit(FatalIOError);
            }
            netStorage_[n] = k;
        }
        else if (roles_[k] == role::fixedPressure)
        {
            if (netStorage_[n] >= 0 || netFixed_[n] >= 0)
            {
                FatalIOErrorInFunction(dict)
                    << "'" << comps_[k].name() << "' shares a node with another"
                    << " compliance, chamber or pressure source. A node can"
                    << " hold only one. Members: " << netMembers(n)
                    << exit(FatalIOError);
            }
            netFixed_[n] = k;
        }
    }

    for (label n = 0; n < nNets_; ++n)
    {
        if (attachments[n] < 2)
        {
            FatalIOErrorInFunction(dict)
                << "Dead end: only " << netMembers(n) << " is attached to a"
                << " node. Every terminal must be connected to something."
                << exit(FatalIOError);
        }

        if (netStorage_[n] < 0 && netFixed_[n] < 0)
        {
            algIndex_[n] = nAlg_++;
        }
    }

    // ---- 3b. Groups of algebraic nets without a pressure reference --------

    groupOfNet_.setSize(nNets_, -1);
    pinned_.setSize(nNets_, false);
    forcedOfComp_.setSize(nComps, -1);

    {
        labelList gp(nNets_);
        std::iota(gp.begin(), gp.end(), 0);

        auto gfind = [&](label a)
        {
            while (gp[a] != a)
            {
                gp[a] = gp[gp[a]];
                a = gp[a];
            }
            return a;
        };

        boolList referenced(nNets_, false);

        forAll(comps_, k)
        {
            if (roles_[k] != role::branch)
            {
                continue;
            }

            const label a = netIn_[k];
            const label b = netOut_[k];

            if (algIndex_[a] >= 0 && algIndex_[b] >= 0)
            {
                gp[gfind(a)] = gfind(b);
            }
        }

        forAll(comps_, k)
        {
            if (roles_[k] != role::branch)
            {
                continue;
            }

            const label a = netIn_[k];
            const label b = netOut_[k];

            if (algIndex_[a] >= 0 && algIndex_[b] < 0)
            {
                referenced[gfind(a)] = true;
            }
            if (algIndex_[b] >= 0 && algIndex_[a] < 0)
            {
                referenced[gfind(b)] = true;
            }
        }

        for (label n = 0; n < nNets_; ++n)
        {
            if (algIndex_[n] >= 0)
            {
                groupOfNet_[n] = gfind(n);
            }
        }

        // A referenced group is solved by its resistors alone. A floating one
        // is held by the single inductor that carries its flow.
        for (label root = 0; root < nNets_; ++root)
        {
            if (algIndex_[root] < 0 || gfind(root) != root || referenced[root])
            {
                continue;
            }

            DynamicList<label> inductors;
            forAll(comps_, k)
            {
                if
                (
                    roles_[k] == role::inertia
                 && (
                        groupOfNet_[netIn_[k]] == root
                     || groupOfNet_[netOut_[k]] == root
                    )
                )
                {
                    inductors.append(k);
                }
            }

            string members;
            forAll(comps_, k)
            {
                if
                (
                    groupOfNet_[netIn_[k]] == root
                 || groupOfNet_[netOut_[k]] == root
                )
                {
                    members +=
                        (members.empty() ? "'" : ", '")
                      + comps_[k].name() + "'";
                }
            }

            if (inductors.size() != 1)
            {
                FatalIOErrorInFunction(dict)
                    << "The pressure of the nodes joining " << members << " is"
                    << " undetermined: they have no resistive path to a"
                    << " compliance, chamber or ground, and "
                    << (inductors.empty() ? "no" : "more than one")
                    << " inductor to carry the flow. Add a compliance or"
                    << " ground, or connect through a resistor."
                    << exit(FatalIOError);
            }

            const label k = inductors[0];
            forcedInertia f;
            f.comp = k;
            f.group = root;
            if (groupOfNet_[netIn_[k]] == root)
            {
                f.pinnedNet = netIn_[k];
                f.otherNet = netOut_[k];
                f.sign = 1;
            }
            else
            {
                f.pinnedNet = netOut_[k];
                f.otherNet = netIn_[k];
                f.sign = -1;
            }

            if (groupOfNet_[f.otherNet] == root)
            {
                FatalIOErrorInFunction(dict)
                    << "Inductor '" << comps_[k].name() << "' has both ends in"
                    << " the same unreferenced group" << exit(FatalIOError);
            }

            forcedOfComp_[k] = forced_.size();
            pinned_[f.pinnedNet] = true;
            forced_.append(f);
        }
    }

    // ---- 4. States ---------------------------------------------------------

    stateOfComp_.setSize(nComps, -1);
    nStates_ = 0;
    forAll(comps_, k)
    {
        if (roles_[k] == role::storage)
        {
            stateOfComp_[k] = nStates_++;
        }
    }
    forAll(comps_, k)
    {
        if (roles_[k] == role::inertia && forcedOfComp_[k] < 0)
        {
            stateOfComp_[k] = nStates_++;
        }
    }

    if (nStates_ == 0)
    {
        FatalIOErrorInFunction(dict)
            << "The network has no state: add a compliance, chamber or inductor"
            << exit(FatalIOError);
    }

    // Branches meeting an algebraic net decide whether one Newton step is exact
    forAll(comps_, k)
    {
        if
        (
            roles_[k] == role::branch
         && (algIndex_[netIn_[k]] >= 0 || algIndex_[netOut_[k]] >= 0)
         && !comps_[k].linear()
        )
        {
            algLinear_ = false;
        }

        if (roles_[k] == role::portSource)
        {
            portIndex_.insert(comps_[k].name(), k);
        }
    }

    // ---- 5. Initial state, columns, solver ---------------------------------

    y_.setSize(nStates_);
    forAll(comps_, k)
    {
        if (roles_[k] == role::storage)
        {
            y_[stateOfComp_[k]] = comps_[k].initialVolume();
        }
        else if (roles_[k] == role::inertia && forcedOfComp_[k] < 0)
        {
            y_[stateOfComp_[k]] = comps_[k].initialFlow();
        }
    }

    pAlg_.setSize(nAlg_, 0);
    Pwork_.setSize(nNets_, 0);

    {
        DynamicList<word> cols;

        wordList probes;
        dict.readIfPresent("probes", probes);
        probeNet_.setSize(probes.size());
        forAll(probes, i)
        {
            const label t = lookup(probes[i]);
            probeNet_[i] = netOfTerminal[t];
        }

        forAll(comps_, k)
        {
            const word& nm = comps_[k].name();
            switch (roles_[k])
            {
                case role::storage:
                    cols.append(nm + ".P");
                    cols.append(nm + ".V");
                    break;
                case role::branch:
                case role::inertia:
                case role::fixedFlow:
                    cols.append(nm + ".Q");
                    break;
                case role::portSource:
                    cols.append(nm + ".Q");
                    cols.append(nm + ".Ppatch");
                    break;
                default:
                    break;
            }
        }
        for (const word& p : probes)
        {
            cols.append(p + ".P");
        }
        columns_ = cols;
    }

    odeDict_ = dict.subOrEmptyDict("solver");

    if (!odeDict_.found("solver"))
    {
        odeDict_.add
        (
            "solver",
            odeDict_.getOrDefault<word>("type", "RKDP45"),
            true
        );
    }
    if (!odeDict_.found("absTol") && !odeDict_.found("absTolComponents"))
    {
        odeDict_.add("absTol", scalar(1e-12));
    }
    if (!odeDict_.found("relTol") && !odeDict_.found("relTolComponents"))
    {
        odeDict_.add("relTol", scalar(1e-8));
    }
    if (!odeDict_.found("maxSteps"))
    {
        odeDict_.add("maxSteps", label(100000));
    }
    dxEst_ = odeDict_.getOrDefault<scalar>("initialStep", 1e-4);

    solver_ = ODESolver::New(*this, odeDict_);
}


Foam::string Foam::circulationNetwork::netMembers(const label net) const
{
    string s;
    forAll(comps_, k)
    {
        if (netIn_[k] == net || netOut_[k] == net)
        {
            if (!s.empty())
            {
                s += ", ";
            }
            s += "'" + comps_[k].name() + "'";
        }
    }
    return s;
}


Foam::scalar Foam::circulationNetwork::branchFlow
(
    const label k,
    const scalar t,
    const scalarField& y,
    const scalarField& P
) const
{
    switch (roles_[k])
    {
        case role::branch:
            return comps_[k].flow(t, P[netIn_[k]] - P[netOut_[k]]);
        case role::inertia:
            return
                forcedOfComp_[k] >= 0
              ? forcedFlow(k, t)
              : y[stateOfComp_[k]];
        case role::fixedFlow:
            return comps_[k].fixedFlow(t);
        default:
            return 0;
    }
}


Foam::scalar Foam::circulationNetwork::forcedFlow
(
    const label k,
    const scalar t
) const
{
    const forcedInertia& f = forced_[forcedOfComp_[k]];

    // Net flow injected into the group by the sources attached to it
    scalar inflow = 0;
    forAll(comps_, s)
    {
        if (roles_[s] == role::fixedFlow)
        {
            const scalar Q = comps_[s].fixedFlow(t);
            if (groupOfNet_[netOut_[s]] == f.group)
            {
                inflow += Q;
            }
            if (groupOfNet_[netIn_[s]] == f.group)
            {
                inflow -= Q;
            }
        }
        else if (roles_[s] == role::portSource)
        {
            if (groupOfNet_[netIn_[s]] == f.group)
            {
                inflow += comps_[s].injectedFlow(t);
            }
        }
    }

    // in -> out through the inductor: it draws the inflow if its pinned end
    // is "in", and delivers it if the pinned end is "out"
    return f.sign*inflow;
}


Foam::scalar Foam::circulationNetwork::forcedFlowRate
(
    const label k,
    const scalar t
) const
{
    const scalar h = 1e-6;
    return (forcedFlow(k, t + h) - forcedFlow(k, t - h))/(2*h);
}


void Foam::circulationNetwork::nodalSystem
(
    const scalar t,
    const scalarField& y,
    scalarField& P,
    const scalarField& pa,
    scalarField& R,
    scalarSquareMatrix* J,
    scalar& flowScale
) const
{
    forAll(algIndex_, n)
    {
        if (algIndex_[n] >= 0)
        {
            P[n] = pa[algIndex_[n]];
        }
    }

    R = Zero;
    if (J)
    {
        *J = Zero;
    }
    flowScale = 0;

    forAll(comps_, k)
    {
        const label ni = netIn_[k];
        const label no = netOut_[k];
        const label ai = algIndex_[ni];
        const label ao = algIndex_[no];

        scalar Q = 0, g = 0;

        switch (roles_[k])
        {
            case role::branch:
            {
                if (ai < 0 && ao < 0)
                {
                    continue;
                }
                const scalar dP = P[ni] - P[no];
                Q = comps_[k].flow(t, dP);
                if (J)
                {
                    g = comps_[k].dFlowdDP(t, dP);
                }
                break;
            }
            case role::inertia:
                Q = branchFlow(k, t, y, P);
                break;
            case role::fixedFlow:
                Q = comps_[k].fixedFlow(t);
                break;
            case role::portSource:
                if (ai >= 0 && !pinned_[ni])
                {
                    const scalar q = comps_[k].injectedFlow(t);
                    R[ai] += q;
                    flowScale = Foam::max(flowScale, Foam::mag(q));
                }
                continue;
            default:
                continue;
        }

        flowScale = Foam::max(flowScale, Foam::mag(Q));

        if (ai >= 0 && !pinned_[ni])
        {
            R[ai] -= Q;
            if (J)
            {
                (*J)(ai, ai) -= g;
                if (ao >= 0)
                {
                    (*J)(ai, ao) += g;
                }
            }
        }
        if (ao >= 0 && !pinned_[no])
        {
            R[ao] += Q;
            if (J)
            {
                (*J)(ao, ao) -= g;
                if (ai >= 0)
                {
                    (*J)(ao, ai) += g;
                }
            }
        }
    }

    // Rows of nets held by a forced inductor: P_in - P_out = L dQ/dt
    for (const forcedInertia& f : forced_)
    {
        const label a = algIndex_[f.pinnedNet];
        const label b = algIndex_[f.otherNet];
        const scalar LdQ =
            comps_[f.comp].inertance()*forcedFlowRate(f.comp, t);

        // pinned end "in":  P_a - P_b = L Q';
        // pinned end "out": P_b - P_a = L Q'
        R[a] = P[f.pinnedNet] - P[f.otherNet] - f.sign*LdQ;
        if (J)
        {
            for (label j = 0; j < nAlg_; ++j)
            {
                (*J)(a, j) = 0;
            }
            (*J)(a, a) = 1;
            if (b >= 0)
            {
                (*J)(a, b) = -1;
            }
        }
    }
}


void Foam::circulationNetwork::netPressures
(
    const scalar t,
    const scalarField& y,
    scalarField& P
) const
{
    forAll(comps_, k)
    {
        if (roles_[k] == role::storage)
        {
            P[netIn_[k]] = comps_[k].pressure(t, y[stateOfComp_[k]]);
        }
        else if (roles_[k] == role::fixedPressure)
        {
            P[netIn_[k]] = comps_[k].fixedPressure(t);
        }
    }

    if (nAlg_ == 0)
    {
        return;
    }

    scalarSquareMatrix J(nAlg_);
    scalarField R(nAlg_);
    scalar flowScale = 0;

    auto singular = []()
    {
        FatalErrorInFunction
            << "The pressure at an algebraic node is undetermined (singular"
            << " nodal system). Every node that holds no volume needs a"
            << " resistive path to a node that does, or to a ground."
            << exit(FatalError);
    };

    if (algLinear_)
    {
        // Only resistors meet at the algebraic nodes: one Newton step is exact
        nodalSystem(t, y, P, pAlg_, R, &J, flowScale);

        scalarField delta(-R);
        if (!solveDense(J, delta))
        {
            singular();
        }

        forAll(pAlg_, a)
        {
            pAlg_[a] += delta[a];
        }
    }
    else
    {
        // A valve (or another nonlinear element) meets an algebraic node: the
        // flows are not linear in the node pressure. Newton with a line
        // search, so that the switching of a hard valve cannot make it cycle,
        // and an error, never a silent wrong answer, if it does not converge.
        auto flowResidual = [&](const scalarField& r)
        {
            scalar m = 0;
            forAll(algIndex_, n)
            {
                if (algIndex_[n] >= 0 && !pinned_[n])
                {
                    m = Foam::max(m, Foam::mag(r[algIndex_[n]]));
                }
            }
            return m;
        };

        nodalSystem(t, y, P, pAlg_, R, &J, flowScale);
        scalar norm = flowResidual(R);

        scalarField trial(nAlg_), Rt(nAlg_);

        for
        (
            label iter = 0;
            iter < 80 && norm > 1e-13*flowScale + 1e-24;
            ++iter
        )
        {
            scalarField delta(-R);
            if (!solveDense(J, delta))
            {
                singular();
            }

            scalar lambda = 1;
            bool improved = false;
            for (label ls = 0; ls < 40; ++ls)
            {
                trial = pAlg_ + lambda*delta;
                scalar fsT = 0;
                nodalSystem(t, y, P, trial, Rt, nullptr, fsT);
                const scalar nt = flowResidual(Rt);

                if (nt < (1 - 1e-4*lambda)*norm || nt <= 1e-13*fsT + 1e-24)
                {
                    improved = true;
                    break;
                }
                lambda *= 0.5;
            }

            if (!improved)
            {
                break;
            }

            pAlg_ = trial;
            nodalSystem(t, y, P, pAlg_, R, &J, flowScale);
            norm = flowResidual(R);
        }

        // The absolute floor is roundoff: a square-root valve that is almost
        // shut turns the rounding of the node pressure (1e-12 relative) into a
        // flow error of order 1e-14 m3/s that no iteration can remove
        if (norm > 1e-8*flowScale + 1e-14)
        {
            FatalErrorInFunction
                << "The pressures of the nodes that hold no volume could not"
                << " be solved at t = " << t << ": the flows into them are out"
                << " of balance by " << norm << " m3/s (largest flow "
                << flowScale << " m3/s)." << nl
                << "    A hard valve (epsilon 0) next to a node without a"
                << " compliance or chamber is the usual cause: give the valve"
                << " a smoothing  epsilon,  or put a compliance at that node."
                << exit(FatalError);
        }
    }

    forAll(algIndex_, n)
    {
        if (algIndex_[n] >= 0)
        {
            P[n] = pAlg_[algIndex_[n]];
        }
    }
}


// * * * * * * * * * * * * * * * ODESystem * * * * * * * * * * * * * * * * * //

void Foam::circulationNetwork::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydx
) const
{
    netPressures(t, y, Pwork_);
    const scalarField& P = Pwork_;

    dydx = Zero;

    forAll(comps_, k)
    {
        const label ni = netIn_[k];
        const label no = netOut_[k];

        switch (roles_[k])
        {
            case role::inertia:
                if (forcedOfComp_[k] < 0)
                {
                    dydx[stateOfComp_[k]] =
                        (P[ni] - P[no])/comps_[k].inertance();
                }
                // its flow also feeds the nets at its ends
                [[fallthrough]];
            case role::branch:
            case role::fixedFlow:
            {
                const scalar Q = branchFlow(k, t, y, P);

                if (netStorage_[ni] >= 0)
                {
                    dydx[stateOfComp_[netStorage_[ni]]] -= Q;
                }
                if (netStorage_[no] >= 0)
                {
                    dydx[stateOfComp_[netStorage_[no]]] += Q;
                }
                break;
            }
            case role::portSource:
                if (netStorage_[ni] >= 0)
                {
                    dydx[stateOfComp_[netStorage_[ni]]] +=
                        comps_[k].injectedFlow(t);
                }
                break;
            default:
                break;
        }
    }
}


void Foam::circulationNetwork::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    const label n = nStates_;

    scalarField f0(n), f1(n);
    derivatives(x, y, f0);

    scalarField yp(y);
    for (label j = 0; j < n; ++j)
    {
        const scalar h = 1e-6*Foam::max(Foam::mag(y[j]), scalar(1e-6));
        yp[j] = y[j] + h;
        derivatives(x, yp, f1);
        yp[j] = y[j];

        for (label i = 0; i < n; ++i)
        {
            dfdy(i, j) = (f1[i] - f0[i])/h;
        }
    }

    const scalar hx = 1e-8*Foam::max(Foam::mag(x), scalar(1));
    derivatives(x + hx, y, f1);
    for (label i = 0; i < n; ++i)
    {
        dfdx[i] = (f1[i] - f0[i])/hx;
    }
}


// * * * * * * * * * * * * * * * Time stepping * * * * * * * * * * * * * * * //

void Foam::circulationNetwork::advance(const scalar t1)
{
    if (t1 <= t_)
    {
        return;
    }

    // Re-raise solver failures with advice. The solver's own error is made
    // catchable for the duration of the call.
    const bool wasThrowing = FatalError.throwing(true);
    try
    {
        solver_->solve(t_, t1, y_, dxEst_);
    }
    catch (Foam::error& e)
    {
        FatalError.throwing(wasThrowing);

        const std::string what = e.message();
        const bool stepCap =
            what.find("Integration steps greater than maximum")
         != std::string::npos;

        FatalErrorInFunction
            << "The circulation ODE solver (" << odeDict_.get<word>("solver")
            << ") failed integrating from t = " << t_ << " to " << t1 << "."
            << nl
            << (stepCap
              ? "    It needed more than maxSteps = "
              + std::to_string(odeDict_.getOrDefault<label>("maxSteps", 0))
              + " steps in this interval. Usual causes:\n"
                "      1. The solution has become unstable at this time step"
                " (a coupled explicit\n"
                "         scheme, or a high coupling order, can diverge): use"
                "  scheme implicit;\n"
                "         order 2;  or a smaller deltaT. The states at failure"
                " are given below:\n"
                "         if they are far from physical this is the cause.\n"
                "      2. The tolerances are very tight for non-smooth valve"
                " switching: raise\n"
                "         maxSteps, or loosen absTol / relTol in the solver"
                " dictionary, or give\n"
                "         the valves a smoothing  epsilon.\n"
              : std::string("    ") + what + "\n").c_str()
            << "    State at failure: " << y_ << nl
            << exit(FatalError);
    }
    FatalError.throwing(wasThrowing);
    t_ = t1;
}


void Foam::circulationNetwork::setState(const scalar t, const scalarField& y)
{
    if (y.size() != nStates_)
    {
        FatalErrorInFunction
            << "State has " << y.size() << " entries, network has " << nStates_
            << exit(FatalError);
    }
    t_ = t;
    y_ = y;
}


// * * * * * * * * * * * * * * * * Results * * * * * * * * * * * * * * * * * //

Foam::stringList Foam::circulationNetwork::columnUnits() const
{
    typedef circulationContext::quantity q;

    stringList units(columns_.size());
    forAll(columns_, i)
    {
        const word& c = columns_[i];
        units[i] =
            c.ends_with(".P") || c.ends_with(".Ppatch")
          ? ctx_.label(q::pressure)
          : c.ends_with(".V") ? ctx_.label(q::volume)
          : ctx_.label(q::flow);
    }
    return units;
}


void Foam::circulationNetwork::outputValues(scalarField& values) const
{
    typedef circulationContext::quantity q;

    values.setSize(columns_.size());

    netPressures(t_, y_, Pwork_);

    const scalar fP = ctx_.factor(q::pressure);
    const scalar fV = ctx_.factor(q::volume);

    label c = 0;
    forAll(comps_, k)
    {
        switch (roles_[k])
        {
            case role::storage:
                values[c++] = Pwork_[netIn_[k]]/fP;
                values[c++] = y_[stateOfComp_[k]]/fV;
                break;
            case role::branch:
            case role::inertia:
            case role::fixedFlow:
                values[c++] = branchFlow(k, t_, y_, Pwork_)/fV;
                break;
            case role::portSource:
                values[c++] = comps_[k].portFlow(t_)/fV;
                values[c++] = comps_[k].patchPressure(Pwork_[netIn_[k]], t_)/fP;
                break;
            default:
                break;
        }
    }

    for (const label n : probeNet_)
    {
        values[c++] = Pwork_[n]/fP;
    }
}


// * * * * * * * * * * * * * * * * Checkpoint * * * * * * * * * * * * * * * * //

Foam::dictionary Foam::circulationNetwork::stateDict() const
{
    dictionary d;
    d.add("time", t_);
    d.add("stepEstimate", dxEst_);

    dictionary s;
    forAll(comps_, k)
    {
        if (stateOfComp_[k] >= 0)
        {
            s.add(word(comps_[k].name()), y_[stateOfComp_[k]]);
        }
    }
    d.add("state", s);
    return d;
}


void Foam::circulationNetwork::readState(const dictionary& d)
{
    const dictionary& s = d.subDict("state");

    scalarField y(y_);
    forAll(comps_, k)
    {
        if (stateOfComp_[k] >= 0)
        {
            y[stateOfComp_[k]] = s.get<scalar>(comps_[k].name());
        }
    }

    setState(d.get<scalar>("time"), y);
    dxEst_ = d.getOrDefault<scalar>("stepEstimate", dxEst_);
}


// * * * * * * * * * * * * * * * * * Ports * * * * * * * * * * * * * * * * * //

Foam::label Foam::circulationNetwork::portLabel(const word& port) const
{
    const auto iter = portIndex_.cfind(port);

    if (!iter.good())
    {
        FatalErrorInFunction
            << "No coupling port '" << port << "'. Ports: "
            << portIndex_.sortedToc() << exit(FatalError);
    }
    return iter.val();
}


void Foam::circulationNetwork::setPortForcing
(
    const word& port,
    const scalar t0,
    const scalar t1,
    const externalPort::forcingCoeffs& c
)
{
    dynamic_cast<externalPort&>(comps_[portLabel(port)]).setForcing(t0, t1, c);
}


Foam::scalar Foam::circulationNetwork::portPressure(const word& port) const
{
    const label k = portLabel(port);

    netPressures(t_, y_, Pwork_);
    return comps_[k].patchPressure(Pwork_[netIn_[k]], t_);
}


bool Foam::circulationNetwork::portInlet(const word& port) const
{
    return dynamic_cast<const externalPort&>(comps_[portLabel(port)]).inlet();
}


const Foam::word& Foam::circulationNetwork::portPatch(const word& port) const
{
    return
        dynamic_cast<const externalPort&>(comps_[portLabel(port)]).patchName();
}


Foam::scalar Foam::circulationNetwork::portFlow
(
    const word& port,
    const scalar t
) const
{
    return comps_[portLabel(port)].portFlow(t);
}


// * * * * * * * * * * * * * * * * * Info * * * * * * * * * * * * * * * * * //

void Foam::circulationNetwork::summary(Ostream& os) const
{
    os  << "circulationNetwork: " << comps_.size() << " elements, "
        << nNets_ << " nodes (" << nAlg_ << " algebraic), "
        << nStates_ << " states" << nl
        << "  units " << ctx_.name() << ", Tcycle " << ctx_.period() << " s"
        << nl
        << "  ODE solver " << odeDict_.get<word>("solver") << nl;

    if (!portIndex_.empty())
    {
        os  << "  coupling ports:";
        for (const word& p : portIndex_.sortedToc())
        {
            os << ' ' << p;
        }
        os  << nl;
    }
}

Foam::label Foam::circulationNetwork::check
(
    Ostream& os,
    const wordList* meshPatches
) const
{
    typedef circulationContext::quantity q;

    label nWarnings = 0;
    label nShown = 0;
    auto warn = [&](const string& msg)
    {
        ++nWarnings;
        if (nShown < 8)
        {
            ++nShown;
            os << "  WARNING: " << msg.c_str() << nl;
        }
    };

    auto roleName = [](const role r) -> const char*
    {
        switch (r)
        {
            case role::composite:     return "composite";
            case role::branch:        return "branch";
            case role::inertia:       return "inertia";
            case role::fixedFlow:     return "flow source";
            case role::storage:       return "storage";
            case role::fixedPressure: return "fixed P";
            case role::portSource:    return "CFD port";
        }
        return "?";
    };

    os  << "Units:    " << ctx_.name() << "  (pressure "
        << ctx_.label(q::pressure).c_str()
        << ", volume " << ctx_.label(q::volume).c_str() << ", time s)" << nl
        << "Period:   Tcycle = " << ctx_.period() << " s" << nl << nl;

    // ---- elements ---------------------------------------------------------
    os  << comps_.size() << " elements:" << nl;
    forAll(comps_, k)
    {
        os  << "  " << comps_[k].name().c_str();
        for (label pad = comps_[k].name().size(); pad < 24; ++pad)
        {
            os << ' ';
        }
        os  << comps_[k].type().c_str();
        for (label pad = comps_[k].type().size(); pad < 18; ++pad)
        {
            os << ' ';
        }
        os  << roleName(roles_[k]) << ",  node n" << netIn_[k];
        if (netOut_[k] != netIn_[k])
        {
            os << " -> n" << netOut_[k];
        }
        os  << nl;
    }

    // ---- nodes ------------------------------------------------------------
    os  << nl << nNets_ << " nodes (" << nAlg_ << " solved algebraically):"
        << nl;
    for (label n = 0; n < nNets_; ++n)
    {
        os  << "  n" << n << "  ";
        if (netStorage_[n] >= 0)
        {
            os  << "stores volume ("
                << comps_[netStorage_[n]].name().c_str() << ")";
        }
        else if (netFixed_[n] >= 0)
        {
            os  << "fixed pressure ("
                << comps_[netFixed_[n]].name().c_str() << ")";
        }
        else
        {
            os  << "algebraic";
            if (pinned_[n])
            {
                os << ", held by an inductor (no pressure reference)";
            }
        }
        os  << ":  " << netMembers(n).c_str() << nl;
    }

    // ---- states -----------------------------------------------------------
    os  << nl << nStates_ << " ODE states:" << nl;
    labelList compOfState(nStates_, -1);
    forAll(stateOfComp_, k)
    {
        if (stateOfComp_[k] >= 0)
        {
            compOfState[stateOfComp_[k]] = k;
        }
    }
    forAll(compOfState, i)
    {
        const label k = compOfState[i];
        const bool volume = roles_[k] == role::storage;
        os  << "  y[" << i << "]  " << comps_[k].name().c_str()
            << (volume ? "  volume  " : "  flow    ")
            << ctx_.fromSI(volume ? q::volume : q::flow, y_[i]) << ' '
            << ctx_.label(volume ? q::volume : q::flow).c_str()
            << "  (initial)" << nl;
    }
    for (const forcedInertia& f : forced_)
    {
        os  << "  (" << comps_[f.comp].name().c_str()
            << " carries the flow of nodes with no"
            << " pressure reference: its flow is not a state, its pressure"
               " drop is L dQ/dt)" << nl;
    }

    // ---- coupling ports ---------------------------------------------------
    if (!portIndex_.empty())
    {
        os  << nl << portIndex_.size() << " coupling port(s):" << nl;
        for (const word& pn : portIndex_.sortedToc())
        {
            const externalPort& ep =
                dynamic_cast<const externalPort&>(comps_[portIndex_[pn]]);
            os  << "  " << pn.c_str() << "  "
                << (
                       ep.inlet()
                     ? "inlet (flow into the CFD domain)"
                     : "outlet (flow leaves the CFD domain)"
                   )
                << ",  patch ";
            if (ep.patchName().empty())
            {
                os  << "(none)" << nl;
                if (!ep.tabulated())
                {
                    warn("port " + pn + " has no \"patch\" entry");
                }
            }
            else
            {
                os  << ep.patchName().c_str();
                if (meshPatches)
                {
                    const bool found = meshPatches->found(ep.patchName());
                    os  << (
                               found
                             ? "  (found in the mesh)"
                             : "  (NOT in the mesh)"
                           )
                        << nl;
                    if (!found)
                    {
                        warn
                        (
                            "patch '" + ep.patchName() + "' of port " + pn
                          + " is not in the mesh"
                        );
                    }
                }
                else
                {
                    os  << nl;
                }
            }
        }
        if (!meshPatches)
        {
            os  << "  (no mesh in this case: patch names not checked)" << nl;
        }
    }

    // ---- parameter magnitudes ---------------------------------------------
    struct range { q quantity; scalar lo, hi; const char* what; };
    static const range ranges[] =
    {
        {q::resistance,        1e3,   1e11,  "resistance"},
        {q::compliance,        1e-12, 1e-5,  "compliance"},
        {q::inertance,         1e2,   1e8,   "inertance"},
        {q::elastance,         1e5,   1e12,  "elastance"},
        {q::volume,            1e-9,  1e-2,  "volume"},
        {q::pressure,          1.0,   1e7,   "pressure"},
        {q::valveCoefficient,  1e-7,  1e-2,  "valve coefficient"},
        {q::inverseVolume,     1e2,   1e9,   "exponential stiffness"},
    };

    os  << nl << "Parameter check (human-scale ranges, SI):" << nl;
    label nChecked = 0;
    forAll(comps_, k)
    {
        DynamicList<circulationComponent::parameter> ps;
        comps_[k].parameters(ps);

        for (const circulationComponent::parameter& p : ps)
        {
            for (const range& r : ranges)
            {
                if (r.quantity != p.quantity)
                {
                    continue;
                }
                ++nChecked;

                const scalar a = Foam::mag(p.value);
                // Zero offsets and unstressed volumes are normal
                if (a == 0)
                {
                    break;
                }

                if (a < r.lo || a > r.hi)
                {
                    std::ostringstream m;
                    m   << comps_[k].name() << "." << p.name << " = "
                        << ctx_.fromSI(p.quantity, p.value) << ' '
                        << ctx_.label(p.quantity) << " (" << p.value
                        << " in SI) is outside the usual "
                        << r.what << " range " << r.lo << " .. " << r.hi
                        << " (SI)";
                    if (ctx_.name() == "SI")
                    {
                        m << ". Written in mmHg and mL? then use  units "
                             "mmHg_mL;";
                    }
                    else
                    {
                        m << ". Is this an SI value in a  units " << ctx_.name()
                          << ";  case?";
                    }
                    warn(m.str());
                }
                break;
            }
        }
    }
    if (nWarnings > nShown)
    {
        os  << "  ... and " << (nWarnings - nShown)
            << " more outside the usual ranges" << nl;
    }
    if (nWarnings == 0)
    {
        os  << "  " << nChecked
            << " values checked, none outside the usual ranges" << nl;
    }

    // ---- solver -----------------------------------------------------------
    {
        std::ostringstream o;
        o.precision(6);
        o   << "ODE solver:  " << odeDict_.get<word>("solver") << ",  absTol "
            << odeDict_.getOrDefault<scalar>("absTol", -1) << ",  relTol "
            << odeDict_.getOrDefault<scalar>("relTol", -1) << ",  maxSteps "
            << odeDict_.getOrDefault<label>("maxSteps", -1);
        os  << nl << o.str().c_str() << nl;
    }

    os  << nl << "Result: " << nStates_ << " states, " << nNets_ << " nodes, "
        << nWarnings << " warning(s)" << nl;

    return nWarnings;
}


// ************************************************************************* //
