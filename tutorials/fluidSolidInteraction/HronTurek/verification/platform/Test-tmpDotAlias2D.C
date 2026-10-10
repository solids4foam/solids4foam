// Two-dimensional variant of the 3dTube reproducer Test-tmpDotAlias.C
// (origin/verification/3dtube-level3, verification/platform/): the same
// tmp<vectorField> & (-tmp<symmTensorField>) expression and loop structure (the
// generated code depends on the context, not on the data), with the inputs of
// a 2-D x-y case: n_z = 0 and D_xz = D_yz = 0, D_zz != 0. The result reuses the
// storage of the first tmp operand while the field loop accesses result and
// operands through __restrict__ pointers. With GCC 11.4 -O3 on xenosim the
// 3-D inputs give a wrong z component only (x, y exact); these 2-D inputs give
// an exact result in all components. Build as buildAliasTest.sh does.
#include "vectorField.H"
#include "symmTensorField.H"
#include "tmp.H"
#include "IOstreams.H"
#include <cmath>

using namespace Foam;

int main()
{
    const label N = 1600;
    vectorField n0(N);
    symmTensorField d0(N);
    forAll(n0, i)
    {
        const scalar a = 0.001*i;
        n0[i] = vector(std::sin(a), std::cos(a), 0.0);
        d0[i] = symmTensor
        (
            1e-4*std::cos(a), 6e-6*std::sin(2*a), 0.0,
            6e-5*std::sin(a), 0.0, 1e-4*std::sin(5*a)
        );
    }

    // 1. As in pimpleFluid::patchViscousForce: tmp & (-tmp)
    tmp<vectorField> tn(new vectorField(n0));
    tmp<symmTensorField> td(new symmTensorField(d0));
    const vectorField viaTmp(tn & (-td));

    // 2. Named operands: no storage reuse
    const vectorField named(n0 & (-d0));

    // 3. By hand
    vectorField byHand(N);
    forAll(byHand, i)
    {
        const vector& n = n0[i];
        const symmTensor d(-d0[i]);
        byHand[i] = vector
        (
            n.x()*d.xx() + n.y()*d.xy() + n.z()*d.xz(),
            n.x()*d.xy() + n.y()*d.yy() + n.z()*d.yz(),
            n.x()*d.xz() + n.y()*d.yz() + n.z()*d.zz()
        );
    }

    vector eC(Zero);
    forAll(byHand, i) { eC = max(eC, cmptMag(viaTmp[i] - byHand[i])); }
    Info<< "per-component max|tmp-reuse - byHand| = " << eC << nl;
    scalar eTmp = 0, eNamed = 0;
    forAll(byHand, i)
    {
        eTmp = max(eTmp, mag(viaTmp[i] - byHand[i]));
        eNamed = max(eNamed, mag(named[i] - byHand[i]));
    }
    Info.stream().precision(17);
    Info<< "max|tmp-reuse - byHand| = " << eTmp
        << "  max|named - byHand| = " << eNamed
        << "  viaTmp[1] = " << viaTmp[1] << "  byHand[1] = " << byHand[1]
        << nl;
    return (eTmp > 1e-12);
}
