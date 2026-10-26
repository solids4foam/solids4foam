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


Application
    Test-tmpAliasProducts

Description
    Checks field products in which a temporary operand is reused as the
    storage of the result.

    OpenFOAM's product operators, e.g. tmp<vectorField> & tmp<tensorField>,
    write the result into the storage of a temporary operand when the types
    allow it (reuseTmp/reuseTmpTmp). The field loops access the result and the
    operands through __restrict__ pointers, so the result then aliases an
    input. For products where each result component depends on several
    components of the reused operand, some GCC -O3 builds give wrong results.
    This was the cause of wrong wall shear stress in the 3dTube FSI case: the
    normal field was an unnamed temporary in "nf() & (-devReff())".

    Each product is evaluated with temporary operands and with named
    operands, which never reuse storage, and the results are compared. The
    application returns a non-zero exit code if any product differs.

    The failure is not in solids4foam but in the OpenFOAM operator templates
    combined with the compiler: it depends on the compiler and on the -O3
    vectorisation, e.g. it is seen with GCC 11 on Ubuntu 22.04 for
    "tmp vector & symmTensor" and not with Clang. A failing line identifies a
    hazardous pattern that must not appear in solids4foam: name the operands
    (e.g. "const vectorField nf(...); nf & x") so that no storage is reused.
    solids4foam should not contain any of the failing patterns.

\*---------------------------------------------------------------------------*/

#include "vectorField.H"
#include "tensorField.H"
#include "symmTensorField.H"
#include "tmp.H"
#include "IOstreams.H"
#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

template<class Type>
scalar maxDiff(const Field<Type>& a, const Field<Type>& b)
{
    scalar e = 0;
    forAll(a, i)
    {
        e = max(e, scalar(mag(a[i] - b[i])));
    }
    return e;
}


template<class Type>
bool report
(
    const char* name,
    const tmp<Field<Type>>& tViaTmp,
    const tmp<Field<Type>>& tRef
)
{
    const scalar e = maxDiff(tViaTmp(), tRef());
    const bool ok = (e < 1e-12);
    Info<< (ok ? "  ok     " : "  FAILED ") << name
        << ": max|tmp - named| = " << e << nl;
    return ok;
}


// * * * * * * * * * * * * * * * * Main  * * * * * * * * * * * * * * * * * * //

int main()
{
    const label N = 1600;

    vectorField v0(N), w0(N);
    tensorField t0(N), u0(N);
    symmTensorField s0(N);

    forAll(v0, i)
    {
        const scalar a = 0.001*i;
        v0[i] = vector(std::sin(a), std::cos(a), 1e-16*std::sin(3*a));
        w0[i] = vector(std::cos(2*a), std::sin(3*a), std::cos(a));
        t0[i] = tensor
        (
            1e-4*std::cos(a), 6e-6*std::sin(2*a), 2e-5*std::cos(3*a),
            6e-5*std::sin(a), 7e-4*std::cos(a), 1e-4*std::sin(5*a),
            3e-5*std::cos(4*a), 2e-4*std::sin(7*a), 5e-4*std::cos(2*a)
        );
        u0[i] = tensor
        (
            std::sin(a), std::cos(2*a), std::sin(3*a),
            std::cos(4*a), std::sin(5*a), std::cos(6*a),
            std::sin(7*a), std::cos(8*a), std::sin(9*a)
        );
        s0[i] = symmTensor
        (
            1e-4*std::cos(a), 6e-6*std::sin(2*a), 2e-5*std::cos(3*a),
            7e-4*std::cos(a), 1e-4*std::sin(5*a), 5e-4*std::cos(2*a)
        );
    }

    // Fresh temporaries, as would arise from e.g. mesh.nf() or devReff()
    #define TV(x) tmp<vectorField>(new vectorField(x))
    #define TT(x) tmp<tensorField>(new tensorField(x))
    #define TS(x) tmp<symmTensorField>(new symmTensorField(x))

    bool ok = true;

    // Check that the named reference is itself correct by comparing the
    // first product with a hand calculation
    {
        vectorField byHand(N);
        forAll(byHand, i)
        {
            const vector& n = v0[i];
            const symmTensor d(-s0[i]);
            byHand[i] = vector
            (
                n.x()*d.xx() + n.y()*d.xy() + n.z()*d.xz(),
                n.x()*d.xy() + n.y()*d.yy() + n.z()*d.yz(),
                n.x()*d.xz() + n.y()*d.yz() + n.z()*d.zz()
            );
        }
        ok &= report
        (
            "named vector & -symmTensor vs hand calculation",
            tmp<vectorField>(new vectorField(v0 & (-s0))),
            tmp<vectorField>(new vectorField(byHand))
        );
    }

    // The original defect: tmp<vector> & (-tmp<symmTensor>)
    ok &=
        report("tmp vector & -tmp symmTensor", TV(v0) & (-TS(s0)), v0 & (-s0));
    ok &= report("tmp vector & tmp tensor", TV(v0) & TT(t0), v0 & t0);
    ok &= report("tmp vector & const tensor", TV(v0) & t0, v0 & t0);
    ok &= report("tmp vector & const symmTensor", TV(v0) & s0, v0 & s0);
    ok &= report("tmp tensor & tmp vector", TT(t0) & TV(v0), t0 & v0);
    ok &= report("const tensor & tmp vector", t0 & TV(v0), t0 & v0);
    ok &= report("tmp tensor & tmp tensor", TT(t0) & TT(u0), t0 & u0);
    ok &= report("tmp tensor & const tensor", TT(t0) & u0, t0 & u0);
    ok &= report("const tensor & tmp tensor", t0 & TT(u0), t0 & u0);
    ok &= report("tmp tensor & const symmTensor", TT(t0) & s0, t0 & s0);
    ok &= report("tmp vector ^ tmp vector", TV(v0) ^ TV(w0), v0 ^ w0);
    ok &= report("tmp vector ^ const vector", TV(v0) ^ w0, v0 ^ w0);
    ok &= report("dev(tmp tensor)", dev(TT(t0)), dev(t0));

    Info<< (ok ? "PASSED" : "FAILED") << nl;

    return ok ? 0 : 1;
}


// ************************************************************************* //
