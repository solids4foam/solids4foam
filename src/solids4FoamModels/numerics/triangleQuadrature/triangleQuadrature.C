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

#include "triangleQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Map<triangleQuadrature::quadratureRule> triangleQuadrature::rules_;

// Initialize and return quadrature rules
void triangleQuadrature::constructRules()
{
    if (rules_.size() != 0)
    {
        FatalErrorInFunction
            << "attempt to re-construct rules when they already exist"
            << exit(FatalError);
    }

    // 1 point quadrature, exact for polynomials up to 1 order
    rules_.insert(1, quadratureRule{
        List<point>{vector(1.0/3.0, 1.0/3.0, 1.0/3.0)},
        List<scalar>{1.0}
    });

    // 3 point quadrature, exact for polynomials up to 2 order
    rules_.insert(3, quadratureRule{
        List<point>
        {
            point(2.0/3.0, 1.0/6.0, 1.0/6.0),
            point(1.0/6.0, 2.0/3.0, 1.0/6.0),
            point(1.0/6.0, 1.0/6.0, 2.0/3.0)
        },
        List<scalar>{1.0/3.0, 1.0/3.0, 1.0/3.0}
    });

    // 4 point quadrature, exact for polynomials up to 3 order
    rules_.insert(4, quadratureRule{
        List<point>
        {
            point(1.0/3.0, 1.0/3.0, 1.0/3.0),
            point(0.6, 0.2, 0.2),
            point(0.2, 0.6, 0.2),
            point(0.2, 0.2, 0.6)
        },
        List<scalar>{-9.0/16.0, 25.0/48.0, 25.0/48.0, 25.0/48.0}
    });

    // 6 point quadrature, exact for polynomials up to 4 order
    rules_.insert(6, quadratureRule{
        List<point>
        {
            point(0.108103018168070, 0.445948490915965, 0.445948490915965),
            point(0.445948490915965, 0.108103018168070, 0.445948490915965),
            point(0.445948490915965, 0.445948490915965, 0.108103018168070),
            point(0.816847572980459, 0.091576213509771, 0.091576213509771),
            point(0.091576213509771, 0.816847572980459, 0.091576213509771),
            point(0.091576213509771, 0.091576213509771, 0.816847572980459)
        },
        List<scalar>
        {
            0.223381589678011,
            0.223381589678011,
            0.223381589678011,
            0.109951743655322,
            0.109951743655322,
            0.109951743655322
        }
    });

    // 7 point quadrature, exact for polynomials up to 5 order
    rules_.insert(7, quadratureRule{
        List<point>
        {
            point(1.0/3.0, 1.0/3.0, 1.0/3.0),
            point(0.059715871789770, 0.470142064105115, 0.470142064105115),
            point(0.470142064105115, 0.059715871789770, 0.470142064105115),
            point(0.470142064105115, 0.470142064105115, 0.059715871789770),
            point(0.797426985353087, 0.101286507323456, 0.101286507323456),
            point(0.101286507323456, 0.797426985353087, 0.101286507323456),
            point(0.101286507323456, 0.101286507323456, 0.797426985353087),
        },
        List<scalar>
        {
             0.225000000000000,
             0.132394152788506,
             0.132394152788506,
             0.132394152788506,
             0.125939180544827,
             0.125939180544827,
             0.125939180544827
        }
    });

    // 12 point quadrature, exact for polynomials up to 6 order
    rules_.insert(12, quadratureRule{
        List<point>
        {
            point(0.501426509658179, 0.249286745170910, 0.249286745170910),
            point(0.249286745170910, 0.501426509658179, 0.249286745170910),
            point(0.249286745170910, 0.249286745170910, 0.501426509658179),
            point(0.873821971016996, 0.063089014491502, 0.063089014491502),
            point(0.063089014491502, 0.873821971016996, 0.063089014491502),
            point(0.063089014491502, 0.063089014491502, 0.873821971016996),
            point(0.053145049844817, 0.310352451033784, 0.636502499121399),
            point(0.053145049844817, 0.636502499121399, 0.310352451033784),
            point(0.310352451033784, 0.053145049844817, 0.636502499121399),
            point(0.636502499121399, 0.053145049844817, 0.310352451033784),
            point(0.310352451033784, 0.636502499121399, 0.053145049844817),
            point(0.636502499121399, 0.310352451033784, 0.053145049844817)
        },
        List<scalar>
        {
            0.116786275726379,
            0.116786275726379,
            0.116786275726379,
            0.050844906370207,
            0.050844906370207,
            0.050844906370207,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374,
            0.082851075618374
        }
    });
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

tmp<Field<point>> triangleQuadrature::barycentricToPoint
(
    const List<point>& localGP
) const
{
    tmp<Field<point>> tglobalGP(new Field<point>(localGP.size()));
    Field<point>& globalGP = tglobalGP.ref();

    forAll(globalGP, pointI)
    {
        globalGP[pointI] =
                localGP[pointI].x() * this->a()
              + localGP[pointI].y() * this->b()
              + localGP[pointI].z() * this->c();
    }

    return tglobalGP;
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Map<triangleQuadrature::quadratureRule>& triangleQuadrature::rules()
{
    if (rules_.size() == 0)
    {
        constructRules();
    }

    return rules_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


triangleQuadrature::triangleQuadrature(const triPoints& pts, const label& n)
:
    triPoints(pts),
    n_(n)
{
    // Check if the requested number of points exists in the rules
    if (!rules().found(n_))
    {
         FatalErrorInFunction
            << "Gaussian quadrature rule for " << n_ << " points not available!"
            << abort(FatalError);
    }

    const quadratureRule& r = rules()[n_];
    weights_ = r.weights;
    gaussPoints_ = barycentricToPoint(r.gaussPoints);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const List<point>& triangleQuadrature::gaussPoints() const
{
    return gaussPoints_;
}


const List<scalar>& triangleQuadrature::weights() const
{
    return weights_;
}


const List<point> triangleQuadrature::gaussPoints()
{
    return gaussPoints_;
}


const List<scalar> triangleQuadrature::weights()
{
     return weights_;
}


} // End namespace Foam

// ************************************************************************* //
