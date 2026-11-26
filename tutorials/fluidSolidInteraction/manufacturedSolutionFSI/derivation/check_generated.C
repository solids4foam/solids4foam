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
    Independent check of the generated header mmsFsiExact.H.

    1. Fluid: finite-difference (4th-order central) evaluation of
         dU/dt + (U.grad)U + grad p - nu lap U   (kinematic, as in solids4foam)
       from the generated U and p only, compared with the generated source f;
       div U; gradU vs FD.
    2. Solid: the neo-Hookean law re-implemented here line by line from
       solids4foam neoHookeanElasticMechanicalConstitutiveLaw.C (evaluate):
         bEbar = J^(-2/3) symm(F F^T);
         sigma = mu/J dev(bEbar) + 0.5 K (J^2-1)/J I
       with F from FD of the generated D; P = J sigma F^-T; then
         rho d2D/dt2 - Div P  (FD)  compared with the generated b0.
    3. Interface: fluid velocity at the deformed interface = solid velocity;
       sigma_f n = sigma_s n with the generated fluid and solid fields.
    4. Magnitudes of every term on a sample grid (how big the sources are
       compared with the physical terms).

    This is a plain C++17 program (no OpenFOAM headers). run_checks.sh
    builds and runs it; to build by hand:

        c++ -O1 -std=c++17 -I ../src \
            -DREFFILE=\"reference_values.txt\" \
            check_generated.C -o check_generated

\*---------------------------------------------------------------------------*/

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "mmsFsiExact.H"

// * * * * * * * * * * * * * * * * Parameters  * * * * * * * * * * * * * * * //

static const double rhoF = mmsFsi::param_rho_f;
static const double muF = mmsFsi::param_mu_f;
static const double nuF = muF/rhoF;
static const double rhoS = mmsFsi::param_rho_s;
static const double Es = mmsFsi::param_E_s;
static const double nus = mmsFsi::param_nu_s;
static const double L = mmsFsi::param_L;
static const double H = mmsFsi::param_H;
static const double T = mmsFsi::param_T;
static int nFail = 0;


// * * * * * * * * * * * * * * * Helper Functions  * * * * * * * * * * * * * //

static void expect(const char* what, double err, double tol)
{
    const bool ok = std::isfinite(err) && err < tol;
    std::printf
    (
        "[%s] %-58s %.2e (tol %.0e)\n", ok ? "PASS" : "FAIL", what, err, tol
    );
    if (!ok)
    {
        nFail++;
    }
}


struct FS
{
    double U[2];
    double p;
    double f[2];
    double G[4];
    double gp[2];
};


static FS fl(double x, double y, double t)
{
    FS r;
    double pp[1];
    mmsFsi::fluid(x, y, t, r.U, pp, r.f, r.G, r.gp);
    r.p = pp[0];
    return r;
}


// 4th-order central first derivative of a scalar functional
template<class F>
double d1(F f, double h)
{
    return (-f(2*h) + 8*f(h) - 8*f(-h) + f(-2*h))/(12*h);
}


template<class F>
double d2(F f, double h)
{
    return
        (-f(2*h) + 16*f(h) - 30*f(0.0) + 16*f(-h) - f(-2*h))/(12*h*h);
}


// * * * * * * * * * * * * Solid law, from solids4foam  * * * * * * * * * * //

struct M3
{
    double a[3][3];
};


static M3 Fof(double X, double Y, double t)
{
    const double h = 1e-5;
    M3 F{};
    for (int i = 0; i < 2; i++)
    {
        F.a[i][0] =
            (i == 0)
          + d1
            (
                [&](double e)
                {
                    double D[2];
                    double U[2];
                    double b[2];
                    double s[6];
                    mmsFsi::solid(X + e, Y, t, D, U, b, s);
                    return D[i];
                },
                h
            );
        F.a[i][1] =
            (i == 1)
          + d1
            (
                [&](double e)
                {
                    double D[2];
                    double U[2];
                    double b[2];
                    double s[6];
                    mmsFsi::solid(X, Y + e, t, D, U, b, s);
                    return D[i];
                },
                h
            );
    }
    F.a[2][2] = 1.0;
    return F;
}


static void neoHookeanSigma(const M3& F, double sig[3][3], double& J)
{
    const double mu = Es/(2*(1 + nus));
    const double lambda = Es*nus/((1 + nus)*(1 - 2*nus));
    const double kappa = lambda + 2.0/3.0*mu;

    // F_zz = 1, plane strain
    J = F.a[0][0]*F.a[1][1] - F.a[0][1]*F.a[1][0];

    double b[3][3] = {};
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                b[i][j] += F.a[i][k]*F.a[j][k];
            }
        }
    }

    const double s = std::pow(J, -2.0/3.0);
    double tr = s*(b[0][0] + b[1][1] + b[2][2]);
    const double sigmaHyd = 0.5*kappa*(J*J - 1.0);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sig[i][j] =
                (mu/J)*(s*b[i][j] - (i == j)*tr/3.0)
              + (i == j)*sigmaHyd/J;
        }
    }
}


static void piola(double X, double Y, double t, double P[2][2])
{
    M3 F = Fof(X, Y, t);
    double sig[3][3];
    double J;
    neoHookeanSigma(F, sig, J);

    // F^-T for the 2x2 block
    const double iT[2][2] =
    {
        { F.a[1][1]/J, -F.a[1][0]/J},
        {-F.a[0][1]/J,  F.a[0][0]/J}
    };

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            P[i][j] = 0;
            for (int k = 0; k < 2; k++)
            {
                P[i][j] += J*sig[i][k]*iT[k][j];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * Main  * * * * * * * * * * * * * * * * //

int main()
{
    std::mt19937 gen(7);
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    double eUD = 0;
    double eF = 0;
    double eDiv = 0;
    double eG = 0;
    double eB = 0;
    double eSig = 0;
    double eKin = 0;
    double eTr = 0;
    double eNrm = 0;
    double fScale = 0;
    double bScale = 0;
    const double h = 1e-4;

    for (int k = 0; k < 200; k++)
    {
        const double t = 0.05 + 0.9*U01(gen);
        const double x = L*(0.02 + 0.96*U01(gen));
        double hh[1];
        mmsFsi::channelHeight(x, t, hh);
        const double y = hh[0]*(0.02 + 0.96*U01(gen));
        FS c = fl(x, y, t);
        double res[2];
        for (int i = 0; i < 2; i++)
        {
            const double dUdt =
                d1([&](double e){ return fl(x, y, t + e).U[i]; }, h);
            const double dUdx =
                d1([&](double e){ return fl(x + e, y, t).U[i]; }, h);
            const double dUdy =
                d1([&](double e){ return fl(x, y + e, t).U[i]; }, h);
            const double lap =
                d2([&](double e){ return fl(x + e, y, t).U[i]; }, 1e-3)
              + d2([&](double e){ return fl(x, y + e, t).U[i]; }, 1e-3);
            const double dp =
                i == 0
              ? d1([&](double e){ return fl(x + e, y, t).p; }, h)
              : d1([&](double e){ return fl(x, y + e, t).p; }, h);
            eG = std::max(eG, std::fabs(dp - c.gp[i]));
            res[i] = dUdt + c.U[0]*dUdx + c.U[1]*dUdy + dp - nuF*lap;
            eF = std::max(eF, std::fabs(res[i] - c.f[i]));
            fScale = std::max(fScale, std::fabs(c.f[i]));
            eG =
                std::max
                (
                    eG,
                    std::fabs(dUdx - c.G[2*i])
                  + std::fabs(dUdy - c.G[2*i + 1])
                );
        }

        // continuity from FD of the generated velocity (not from gradU)
        eDiv =
            std::max
            (
                eDiv,
                std::fabs
                (
                    d1([&](double e){ return fl(x + e, y, t).U[0]; }, h)
                  + d1([&](double e){ return fl(x, y + e, t).U[1]; }, h)
                )
            );

        // solid
        const double X = L*(0.02 + 0.96*U01(gen));
        const double Y = H + T*(0.02 + 0.96*U01(gen));
        double D[2];
        double Us[2];
        double b0[2];
        double sg[6];
        mmsFsi::solid(X, Y, t, D, Us, b0, sg);
        double sig[3][3];
        double J;
        neoHookeanSigma(Fof(X, Y, t), sig, J);
        const double sgRef[6] =
        {
            sig[0][0], sig[0][1], sig[0][2], sig[1][1], sig[1][2], sig[2][2]
        };
        for (int q = 0; q < 6; q++)
        {
            eSig = std::max(eSig, std::fabs(sg[q] - sgRef[q]));
        }
        for (int i = 0; i < 2; i++)
        {
            // U = dD/dt and acceleration = d2D/dt2, both from FD of D only
            const double Dt =
                d1
                (
                    [&](double e)
                    {
                        double D_[2];
                        double U_[2];
                        double b_[2];
                        double s_[6];
                        mmsFsi::solid(X, Y, t + e, D_, U_, b_, s_);
                        return D_[i];
                    },
                    h
                );
            eUD = std::max(eUD, std::fabs(Dt - Us[i]));
            const double acc =
                d2
                (
                    [&](double e)
                    {
                        double D_[2];
                        double U_[2];
                        double b_[2];
                        double s_[6];
                        mmsFsi::solid(X, Y, t + e, D_, U_, b_, s_);
                        return D_[i];
                    },
                    1e-3
                );
            const double divP =
                d1
                (
                    [&](double e)
                    {
                        double P[2][2];
                        piola(X + e, Y, t, P);
                        return P[i][0];
                    },
                    1e-3
                )
              + d1
                (
                    [&](double e)
                    {
                        double P[2][2];
                        piola(X, Y + e, t, P);
                        return P[i][1];
                    },
                    1e-3
                );
            eB = std::max(eB, std::fabs(rhoS*acc - divP - b0[i]));
            bScale = std::max(bScale, std::fabs(b0[i]));
        }

        // interface
        const double XG = L*(0.01 + 0.98*U01(gen));
        double xg[2];
        double tr[2];
        double n[2];
        mmsFsi::interface(XG, t, xg, tr, n);
        mmsFsi::solid(XG, H, t, D, Us, b0, sg);
        FS ci = fl(xg[0], xg[1], t);
        eKin =
            std::max
            (
                eKin,
                std::fabs(ci.U[0] - Us[0]) + std::fabs(ci.U[1] - Us[1])
            );
        for (int i = 0; i < 2; i++)
        {
            const double xt =
                d1
                (
                    [&](double e)
                    {
                        double a_[2];
                        double b_[2];
                        double c_[2];
                        mmsFsi::interface(XG, t + e, a_, b_, c_);
                        return a_[i];
                    },
                    h
                );
            eKin = std::max(eKin, std::fabs(ci.U[i] - xt));
        }
        {
            const double tx =
                d1
                (
                    [&](double e)
                    {
                        double a_[2];
                        double b_[2];
                        double c_[2];
                        mmsFsi::interface(XG + e, t, a_, b_, c_);
                        return a_[0];
                    },
                    h
                );
            const double ty =
                d1
                (
                    [&](double e)
                    {
                        double a_[2];
                        double b_[2];
                        double c_[2];
                        mmsFsi::interface(XG + e, t, a_, b_, c_);
                        return a_[1];
                    },
                    h
                );
            const double m = std::hypot(tx, ty);
            eNrm =
                std::max
                (
                    eNrm,
                    std::fabs(n[0] + ty/m) + std::fabs(n[1] - tx/m)
                );
            double hh2[1];
            mmsFsi::channelHeight(xg[0], t, hh2);
            eNrm = std::max(eNrm, std::fabs(hh2[0] - xg[1]));
        }

        // fluid traction (dynamic): rho*( -p n + nu (G + G^T) n )
        const double G[2][2] = {{ci.G[0], ci.G[1]}, {ci.G[2], ci.G[3]}};
        const double sS[2][2] = {{sg[0], sg[1]}, {sg[1], sg[3]}};
        for (int i = 0; i < 2; i++)
        {
            double tf = -ci.p*n[i];
            double ts = 0;
            for (int j = 0; j < 2; j++)
            {
                tf += nuF*(G[i][j] + G[j][i])*n[j];
                ts += sS[i][j]*n[j];
            }
            tf *= rhoF;
            eTr =
                std::max
                (
                    eTr,
                    std::fabs(tf - ts) + std::fabs(ts - tr[i])
                );
        }
    }

    std::printf
    (
        "FD-based checks (200 random points; "
        "FD truncation ~1e-9..1e-7):\n"
    );
    expect("fluid: FD residual of NS from U,p  vs  generated f", eF, 1e-7);
    expect("fluid: generated gradU, gradp vs FD", eG, 1e-7);
    expect("fluid: div U from FD of U", eDiv, 1e-8);
    expect("solid: generated sigma vs re-implemented s4f law", eSig, 1e-7);
    expect("solid: U vs FD dD/dt", eUD, 1e-8);
    expect
    (
        "solid: rho D_tt - Div P (FD, s4f law) vs generated b0", eB, 1e-5
    );
    expect("interface: v_f(x_G) = U_s(X,H) = d x_G/dt", eKin, 1e-8);
    expect("interface: normal from FD tangent; h(x_G) = y_G", eNrm, 1e-8);
    expect
    (
        "interface: sigma_f n = sigma_s n = generated traction", eTr, 1e-11
    );
    std::printf
    (
        "fluid momentum source  : max |FD residual - f| = %.2e   "
        "(max |f| = %.2e)\n",
        eF,
        fScale
    );
    std::printf("fluid gradU vs FD      : %.2e\n", eG);
    std::printf("fluid div U            : %.2e\n", eDiv);
    std::printf("solid sigma vs s4f law : %.2e\n", eSig);
    std::printf
    (
        "solid source           : max |FD residual - b0| = %.2e   "
        "(max |b0| = %.2e)\n",
        eB,
        bScale
    );
    std::printf("interface kinematics   : %.2e\n", eKin);
    std::printf("interface traction     : %.2e\n", eTr);


    // * * * * * * * * * * * Magnitudes on a grid, several times * * * * * * //

    const int N = 60;
    double mU = 0;
    double mP = 0;
    double mF = 0;
    double mA = 0;
    double mC = 0;
    double mGP = 0;
    double mV = 0;
    double mD = 0;
    double mB = 0;
    double mAcc = 0;
    double mDivP = 0;
    double mStr = 0;
    double minJ = 1e9;
    double maxJ = 0;
    double mRot = 0;
    for (double t : {0.15, 0.3, 0.45, 0.6, 0.75, 0.9})
    {
        for (int i = 1; i < N; i++)
        {
            for (int j = 1; j < N; j++)
            {
                const double x = L*i/N;
                double hh[1];
                mmsFsi::channelHeight(x, t, hh);
                const double y = hh[0]*j/N;
                FS c = fl(x, y, t);
                mU = std::max(mU, std::hypot(c.U[0], c.U[1]));
                mP = std::max(mP, std::fabs(c.p));
                mF = std::max(mF, std::hypot(c.f[0], c.f[1]));
                double a[2];
                double cv[2];
                double gp[2];
                double vs[2];
                for (int q = 0; q < 2; q++)
                {
                    a[q] =
                        d1
                        (
                            [&](double e){ return fl(x, y, t + e).U[q]; },
                            h
                        );
                    cv[q] = c.U[0]*c.G[2*q] + c.U[1]*c.G[2*q + 1];
                    gp[q] =
                        q == 0
                      ? d1([&](double e){ return fl(x + e, y, t).p; }, h)
                      : d1([&](double e){ return fl(x, y + e, t).p; }, h);
                    vs[q] =
                        nuF
                       *(
                            d2
                            (
                                [&](double e)
                                {
                                    return fl(x + e, y, t).U[q];
                                },
                                1e-3
                            )
                          + d2
                            (
                                [&](double e)
                                {
                                    return fl(x, y + e, t).U[q];
                                },
                                1e-3
                            )
                        );
                }
                mA = std::max(mA, std::hypot(a[0], a[1]));
                mC = std::max(mC, std::hypot(cv[0], cv[1]));
                mGP = std::max(mGP, std::hypot(gp[0], gp[1]));
                mV = std::max(mV, std::hypot(vs[0], vs[1]));

                const double X = L*i/N;
                const double Y = H + T*j/N;
                double D[2];
                double Us[2];
                double b0[2];
                double sg[6];
                mmsFsi::solid(X, Y, t, D, Us, b0, sg);
                mD = std::max(mD, std::hypot(D[0], D[1]));
                mB = std::max(mB, std::hypot(b0[0], b0[1]));
                M3 F = Fof(X, Y, t);
                const double J =
                    F.a[0][0]*F.a[1][1] - F.a[0][1]*F.a[1][0];
                minJ = std::min(minJ, J);
                maxJ = std::max(maxJ, J);

                // Green-Lagrange strain magnitude and rotation angle
                // (polar decomposition, 2-D)
                const double C00 =
                    F.a[0][0]*F.a[0][0] + F.a[1][0]*F.a[1][0];
                const double C11 =
                    F.a[0][1]*F.a[0][1] + F.a[1][1]*F.a[1][1];
                const double C01 =
                    F.a[0][0]*F.a[0][1] + F.a[1][0]*F.a[1][1];
                const double E00 = 0.5*(C00 - 1);
                const double E11 = 0.5*(C11 - 1);
                const double E01 = 0.5*C01;
                mStr =
                    std::max
                    (
                        mStr,
                        std::sqrt(E00*E00 + E11*E11 + 2*E01*E01)
                    );
                mRot =
                    std::max
                    (
                        mRot,
                        std::fabs
                        (
                            std::atan2
                            (
                                F.a[1][0] - F.a[0][1],
                                F.a[0][0] + F.a[1][1]
                            )
                        )
                    );
                double acc[2];
                for (int q = 0; q < 2; q++)
                {
                    acc[q] =
                        d1
                        (
                            [&](double e)
                            {
                                double D_[2];
                                double U_[2];
                                double b_[2];
                                double s_[6];
                                mmsFsi::solid
                                (
                                    X, Y, t + e, D_, U_, b_, s_
                                );
                                return U_[q];
                            },
                            h
                        );
                }
                mAcc = std::max(mAcc, rhoS*std::hypot(acc[0], acc[1]));
                mDivP =
                    std::max
                    (
                        mDivP,
                        std::hypot
                        (
                            rhoS*acc[0] - b0[0],
                            rhoS*acc[1] - b0[1]
                        )
                    );
            }
        }
    }
    std::printf("\nmagnitudes (max over sample grid, t in 0.15..0.9):\n");
    std::printf
    (
        " fluid |U| %.3g  |p/rho| %.3g  |f| %.3g | dU/dt %.3g  conv %.3g  "
        "grad p %.3g  visc %.3g\n",
        mU, mP, mF, mA, mC, mGP, mV
    );
    std::printf
    (
        " solid |D| %.3g  |E_GL| %.3g  rot %.3g rad  J in [%.3f, %.3f]  "
        "|b0| %.3g | rho*acc %.3g  |Div P| %.3g\n",
        mD, mStr, mRot, minJ, maxJ, mB, mAcc, mDivP
    );


    // * * * * * * Direct comparison with 40-digit symengine values * * * * * //

    {
        std::ifstream in(REFFILE);
        std::string line;
        double eRef = 0;
        int nRef = 0;
        while (std::getline(in, line))
        {
            std::istringstream ss(line);
            double x;
            double y;
            double t;
            double Y;
            ss >> x >> y >> t >> Y;
            std::vector<double> ref;
            double q;
            while (ss >> q)
            {
                ref.push_back(q);
            }
            double U[2];
            double p[1];
            double f[2];
            double G[4];
            double gp[2];
            double D[2];
            double Us[2];
            double b0[2];
            double sg[6];
            mmsFsi::fluid(x, y, t, U, p, f, G, gp);
            mmsFsi::solid(x, Y, t, D, Us, b0, sg);
            const double got[11] =
            {
                U[0], U[1], p[0], f[0], f[1], D[0], D[1],
                Us[0], Us[1], b0[0], b0[1]
            };
            for (int k = 0; k < 11 && k < int(ref.size()); k++)
            {
                eRef =
                    std::max
                    (
                        eRef,
                        std::fabs(got[k] - ref[k])/(1.0 + std::fabs(ref[k]))
                    );
            }
            nRef++;
        }
        expect
        (
            nRef
          ? "generated double vs 40-digit symbolic values (rel.)"
          : "reference file missing",
            nRef ? eRef : 1.0,
            1e-12
        );
    }
    std::printf
    (
        "%s\n",
        nFail
      ? "SOME CHECKS FAILED"
      : "all generated-code checks passed"
    );
    return nFail ? 1 : 0;
}


// ************************************************************************* //
