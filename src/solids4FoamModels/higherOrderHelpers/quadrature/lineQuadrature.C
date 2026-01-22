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

\*---------------------------------------------------------------------------*/

#include "lineQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Map<lineQuadrature::quadratureRule> lineQuadrature::rules_;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// Initialize and return quadrature rules
void lineQuadrature::constructRules()
{
    if (rules_.size() != 0)
    {
        FatalErrorInFunction
            << "attempt to re-construct rules when they already exist"
            << exit(FatalError);
    }

    // 1 point quadrature, exact for polynomials up to 1 order
    rules_.insert
    (
        1,
        quadratureRule
        {
            List<scalar>{0.0},
            List<scalar>{2.0}
        }
    );

    // 2 point quadrature, exact for polynomials up to 3 order
    rules_.insert
    (
        3,
        quadratureRule
        {
            List<scalar>{1.0/sqrt(3.0), -1.0/sqrt(3.0)},
            List<scalar>{1.0, 1.0}
        }
    );

    // 3 point quadrature, exact for polynomials up to 5 order
    rules_.insert
    (
        5,
        quadratureRule
        {
            List<scalar>
            {
                0.0,
                sqrt(3.0/5.0),
               -sqrt(3.0/5.0)
            },
            List<scalar>{8.0/9.0, 5.0/9.0, 5.0/9.0}
        }
    );

    // 4 point quadrature, exact for polynomials up to 7 order
    rules_.insert
    (
        7,
        quadratureRule
        {
            List<scalar>
            {
                sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0)),
                -sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0)),
                sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0)),
                -sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0))
            },
            List<scalar>
            {
                (18.0+sqrt(30.0))/36.0,
                (18.0+sqrt(30.0))/36.0,
                (18.0-sqrt(30.0))/36.0,
                (18.0-sqrt(30.0))/36.0
            }
        }
    );

    // 5 point quadrature, exact for polynomials up to 9 order
    rules_.insert
    (
        9,
        quadratureRule
        {
            List<scalar>
            {
                0,
                (1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0)),
                -(1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0)),
                (1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0)),
                -(1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0)),
            },
            List<scalar>
            {
                128.0/225.0,
                (322.0+13.0*sqrt(70.0))/900.0,
                (322.0+13.0*sqrt(70.0))/900.0,
                (322.0-13.0*sqrt(70.0))/900.0,
                (322.0-13.0*sqrt(70.0))/900.0
            }
        }
    );

    // 6 point quadrature, exact for polynomials up to 11 order
    rules_.insert
    (
        11,
        quadratureRule
        {
            List<scalar>
            {
                 0.932469514203152,
                -0.932469514203152,
                 0.661209386466265,
                -0.661209386466265,
                 0.238619186083197,
                -0.238619186083197
            },
            List<scalar>
            {
                0.171324492379170,
                0.171324492379170,
                0.360761573048139,
                0.360761573048139,
                0.467913934572691,
                0.467913934572691
            }
        }
    );

    // 7 point quadrature, exact for polynomials up to 13 order
    rules_.insert
    (
        13,
        quadratureRule
        {
            List<scalar>
            {
                0.0,
                0.949107912342759,
               -0.949107912342759,
                0.741531185599394,
               -0.741531185599394,
                0.405845151377397,
               -0.405845151377397,
            },
            List<scalar>
            {
                0.417959183673469,
                0.129484966168870,
                0.129484966168870,
                0.279705391489277,
                0.279705391489277,
                0.381830050505119,
                0.381830050505119
            }
        }
    );

    // 8 point quadrature, exact for polynomials up to 15 order
    rules_.insert
    (
        15,
        quadratureRule
        {
            List<scalar>
            {
                0.960289856497536,
               -0.960289856497536,
                0.796666477413627,
               -0.796666477413627,
                0.525532409916329,
               -0.525532409916329,
                0.183434642495650,
               -0.183434642495650
            },
            List<scalar>
            {
                0.101228536378362,
                0.101228536378362,
                0.222381034453374,
                0.222381034453374,
                0.313706645877887,
                0.313706645877887,
                0.362683783378362,
                0.362683783378362
            }
        }
    );

    // 9 point quadrature, exact for polynomials up to 17 order
    rules_.insert
    (
        17,
        quadratureRule
        {
            List<scalar>
            {
                0.968160239507626,
               -0.968160239507626,
                0.836031107326636,
               -0.836031107326636,
                0.613371432700590,
               -0.613371432700590,
                0.324253423403809,
               -0.324253423403809,
                0.0
            },
            List<scalar>
            {
                0.081274388361574,
                0.081274388361574,
                0.180648160694857,
                0.180648160694857,
                0.260610696402935,
                0.260610696402935,
                0.312347077040003,
                0.312347077040003,
                0.330239355001260
            }
        }
    );

    // 10 point quadrature, exact for polynomials up to 19 order
    rules_.insert
    (
        19,
        quadratureRule
        {
            List<scalar>
            {
                0.973906528517172,
               -0.973906528517172,
                0.865063366688985,
               -0.865063366688985,
                0.679409568299024,
               -0.679409568299024,
                0.433395394129247,
               -0.433395394129247,
                0.148874338981631,
               -0.148874338981631,
            },
            List<scalar>
            {
                0.066671344308688,
                0.066671344308688,
                0.149451349150581,
                0.149451349150581,
                0.219086362515982,
                0.219086362515982,
                0.269266719309996,
                0.269266719309996,
                0.295524224714753,
                0.295524224714753
            }
        }
    );
}


const Map<lineQuadrature::quadratureRule>& lineQuadrature::rules()
{
    if (rules_.size() == 0)
    {
        constructRules();
    }

    return rules_;
}


label lineQuadrature::chooseOrder(const label requestedOrder)
{
    label chosenOrder = -1;

    // If requested order does not exist in the rules, take the first higher
    if (rules().found(requestedOrder))
    {
        chosenOrder = requestedOrder;
    }
    else
    {
        for (label o = requestedOrder + 1; o <= maxSupportedOrder; ++o)
        {
            if (rules().found(o))
            {
                chosenOrder = o;
                break;
            }
        }
    }

    if (chosenOrder == -1)
    {
        FatalErrorInFunction
            << "Quadrature for " << requestedOrder << " order not implemented. "
            << "The highest order of accuracy is " << maxSupportedOrder
            << abort(FatalError);
    }

    return chosenOrder;
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

tmp<Field<point>> lineQuadrature::parametricToPoint
(
    const List<scalar>& x
) const
{
    tmp<Field<point>> tglobalPts(new Field<point>(x.size()));
    Field<point>& globalPts = tglobalPts.ref();

    forAll(globalPts, pointI)
    {
        // Map from [-1,1] to [0,1]
        const scalar t = 0.5*x[pointI] + 0.5;
        globalPts[pointI] = line_.start() + t * line_.vec();
    }

    return tglobalPts;
}


tmp<Field<scalar>> lineQuadrature::normaliseWeights
(
    const List<scalar>& w
) const
{
    tmp<Field<scalar>> tNormWeights(new Field<scalar>(w.size()));
    Field<scalar>& normWeights = tNormWeights.ref();

    // Weights sum is 2 because interval is [-1,1]. When working on
    // interval [0,1] weights should be rescaled.
    forAll(normWeights, wI)
    {
        normWeights[wI] = w[wI]*0.5;
    }

    return tNormWeights;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


lineQuadrature::lineQuadrature
(
    const point& start,
    const point& end,
    const label& order
)
:
    quadrature(order),
    line_(start, end)
{
    const label chosen = chooseOrder(order);
    const auto& rule = rules()[chosen];

    weights_ = rule.weights;
    points_ = parametricToPoint(rule.points);
}

lineQuadrature::lineQuadrature
(
    const line<point, const point&>& l,
    const label& order
)
:
    quadrature(order),
    line_(l.start(), l.end())
{
    const label chosen = chooseOrder(order);
    const auto& rule = rules()[chosen];

    weights_ = rule.weights;
    points_ = parametricToPoint(rule.points);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const List<point>& lineQuadrature::points() const
{
    return points_;
}


const List<scalar>& lineQuadrature::weights() const
{
    return weights_;
}


label lineQuadrature::nPoints() const
{
    return points_.size();
}


label lineQuadrature::nPoints(label order)
{
    const auto& rule = rules()[chooseOrder(order)];
    return rule.points.size();
}

} // End namespace Foam

// ************************************************************************* //
