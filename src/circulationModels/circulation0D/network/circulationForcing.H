/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Description
    The flow forcing a coupling scheme hands to the circulation network over
    one step [t0, t1]: a polynomial of at most third degree in
    alpha = (t - t0)/(t1 - t0),

        Q(alpha) = c0 + c1 alpha + c2 alpha^2 + c3 alpha^3,

    built through up to four (alpha, Q) samples. Samples at alpha < 0 are
    earlier flows (extrapolation, explicit coupling); alpha = 1 is the flow at
    the end of the step (implicit coupling).

\*---------------------------------------------------------------------------*/

#ifndef circulationForcing_H
#define circulationForcing_H

#include "FixedList.H"
#include "UList.H"
#include "scalarList.H"
#include "error.H"

namespace Foam
{
namespace circulationForcing
{

typedef FixedList<scalar, 4> coeffs;

//- Coefficients of the polynomial through n = 1..4 points
inline coeffs polynomialThrough
(
    const UList<scalar>& alpha,
    const UList<scalar>& values
)
{
    const label n = alpha.size();

    if (n < 1 || n > 4 || values.size() != n)
    {
        FatalErrorInFunction
            << "Need one to four forcing points with one value each, given "
            << n << " abscissae and " << values.size() << " values"
            << exit(FatalError);
    }

    // Newton divided differences, in place
    scalarList d(values);
    for (label j = 1; j < n; ++j)
    {
        for (label i = n - 1; i >= j; --i)
        {
            d[i] = (d[i] - d[i-1])/(alpha[i] - alpha[i-j]);
        }
    }

    // Expand into ascending powers: c <- c*(a - a_k) + d_k
    coeffs c(0.0);
    c[0] = d[n-1];
    for (label k = n - 2; k >= 0; --k)
    {
        for (label i = n - 1; i > 0; --i)
        {
            c[i] = c[i-1] - alpha[k]*c[i];
        }
        c[0] = d[k] - alpha[k]*c[0];
    }

    return c;
}

//- Horner evaluation
inline scalar evaluate(const coeffs& c, const scalar alpha)
{
    return c[0] + alpha*(c[1] + alpha*(c[2] + alpha*c[3]));
}

} // End namespace circulationForcing
} // End namespace Foam

#endif
