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

#include "triQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Map<triQuadrature::quadratureRule> triQuadrature::rules_;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// Initialize and return quadrature rules
void triQuadrature::constructRules()
{
    if (rules_.size() != 0)
    {
        FatalErrorInFunction
            << "attempt to re-construct rules when they already exist"
            << exit(FatalError);
    }

    // 1 barycentric2D quadrature, exact for polynomials up to 1 order
    rules_.insert
    (
        1,
        quadratureRule
        {
            List<barycentric2D>{barycentric2D(1.0/3.0, 1.0/3.0, 1.0/3.0)},
            List<scalar>{1.0}
        }
    );

    // 3 barycentric2D quadrature, exact for polynomials up to 2 order
    rules_.insert
    (
        2,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(2.0/3.0, 1.0/6.0, 1.0/6.0),
                barycentric2D(1.0/6.0, 2.0/3.0, 1.0/6.0),
                barycentric2D(1.0/6.0, 1.0/6.0, 2.0/3.0)
            },
            List<scalar>{1.0/3.0, 1.0/3.0, 1.0/3.0}
        }
    );

    // 4 barycentric2D quadrature, exact for polynomials up to 3 order
    rules_.insert
    (
        3,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(1.0/3.0, 1.0/3.0, 1.0/3.0),
                barycentric2D(0.6, 0.2, 0.2),
                barycentric2D(0.2, 0.6, 0.2),
                barycentric2D(0.2, 0.2, 0.6)
            },
            List<scalar>{-9.0/16.0, 25.0/48.0, 25.0/48.0, 25.0/48.0}
        }
    );

    // 6 barycentric2D quadrature, exact for polynomials up to 4 order
    rules_.insert
    (
        4,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(0.108103018168070, 0.445948490915965, 0.445948490915965),
                barycentric2D(0.445948490915965, 0.108103018168070, 0.445948490915965),
                barycentric2D(0.445948490915965, 0.445948490915965, 0.108103018168070),
                barycentric2D(0.816847572980459, 0.091576213509771, 0.091576213509771),
                barycentric2D(0.091576213509771, 0.816847572980459, 0.091576213509771),
                barycentric2D(0.091576213509771, 0.091576213509771, 0.816847572980459)
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
        }
    );

    // 7 barycentric2D quadrature, exact for polynomials up to 5 order
    rules_.insert
    (
        5,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(1.0/3.0, 1.0/3.0, 1.0/3.0),
                barycentric2D(0.059715871789770, 0.470142064105115, 0.470142064105115),
                barycentric2D(0.470142064105115, 0.059715871789770, 0.470142064105115),
                barycentric2D(0.470142064105115, 0.470142064105115, 0.059715871789770),
                barycentric2D(0.797426985353087, 0.101286507323456, 0.101286507323456),
                barycentric2D(0.101286507323456, 0.797426985353087, 0.101286507323456),
                barycentric2D(0.101286507323456, 0.101286507323456, 0.797426985353087),
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
        }
    );

    // 12 barycentric2D quadrature, exact for polynomials up to 6 order
    rules_.insert
    (
        6,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(0.501426509658179, 0.249286745170910, 0.249286745170910),
                barycentric2D(0.249286745170910, 0.501426509658179, 0.249286745170910),
                barycentric2D(0.249286745170910, 0.249286745170910, 0.501426509658179),
                barycentric2D(0.873821971016996, 0.063089014491502, 0.063089014491502),
                barycentric2D(0.063089014491502, 0.873821971016996, 0.063089014491502),
                barycentric2D(0.063089014491502, 0.063089014491502, 0.873821971016996),
                barycentric2D(0.053145049844817, 0.310352451033784, 0.636502499121399),
                barycentric2D(0.053145049844817, 0.636502499121399, 0.310352451033784),
                barycentric2D(0.310352451033784, 0.053145049844817, 0.636502499121399),
                barycentric2D(0.636502499121399, 0.053145049844817, 0.310352451033784),
                barycentric2D(0.310352451033784, 0.636502499121399, 0.053145049844817),
                barycentric2D(0.636502499121399, 0.310352451033784, 0.053145049844817)
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
        }
    );

    // 13 barycentric2D quadrature, exact for polynomials up to 7 order
    rules_.insert
    (
        7,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(0.333333333333333, 0.333333333333333, 0.333333333333333),
                barycentric2D(0.479308067841920, 0.260345966079040, 0.260345966079040),
                barycentric2D(0.260345966079040, 0.479308067841920, 0.260345966079040),
                barycentric2D(0.260345966079040, 0.260345966079040, 0.479308067841920),
                barycentric2D(0.869739794195568, 0.065130102902216, 0.065130102902216),
                barycentric2D(0.065130102902216, 0.869739794195568, 0.065130102902216),
                barycentric2D(0.065130102902216, 0.065130102902216, 0.869739794195568),
                barycentric2D(0.048690315425316, 0.312865496004874, 0.638444188569810),
                barycentric2D(0.048690315425316, 0.638444188569810, 0.312865496004874),
                barycentric2D(0.312865496004874, 0.048690315425316, 0.638444188569810),
                barycentric2D(0.312865496004874, 0.638444188569810, 0.048690315425316),
                barycentric2D(0.638444188569810, 0.048690315425316, 0.312865496004874),
                barycentric2D(0.638444188569810, 0.312865496004874, 0.048690315425316)
            },
            List<scalar>
            {
               -0.149570044467682,
                0.175615257433208,
                0.175615257433208,
                0.175615257433208,
                0.053347235608838,
                0.053347235608838,
                0.053347235608838,
                0.077113760890257,
                0.077113760890257,
                0.077113760890257,
                0.077113760890257,
                0.077113760890257,
                0.077113760890257
            }
        }
    );

    // 16 barycentric2D quadrature, exact for polynomials up to 8 order
    rules_.insert
    (
        8,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(0.333333333333333, 0.333333333333333, 0.333333333333333),
                barycentric2D(0.081414823414554, 0.459292588292723, 0.459292588292723),
                barycentric2D(0.459292588292723, 0.081414823414554, 0.459292588292723),
                barycentric2D(0.459292588292723, 0.459292588292723, 0.081414823414554),
                barycentric2D(0.658861384496480, 0.170569307751760, 0.170569307751760),
                barycentric2D(0.170569307751760, 0.658861384496480, 0.170569307751760),
                barycentric2D(0.170569307751760, 0.170569307751760, 0.658861384496480),
                barycentric2D(0.898905543365938, 0.050547228317031, 0.050547228317031),
                barycentric2D(0.050547228317031, 0.898905543365938, 0.050547228317031),
                barycentric2D(0.050547228317031, 0.050547228317031, 0.898905543365938),
                barycentric2D(0.008394777409958, 0.263112829634638, 0.728492392955404),
                barycentric2D(0.008394777409958, 0.728492392955404, 0.263112829634638),
                barycentric2D(0.263112829634638, 0.008394777409958, 0.728492392955404),
                barycentric2D(0.263112829634638, 0.728492392955404, 0.008394777409958),
                barycentric2D(0.728492392955404, 0.008394777409958, 0.263112829634638),
                barycentric2D(0.728492392955404, 0.263112829634638, 0.008394777409958)
            },
            List<scalar>
            {
                0.144315607677787,
                0.095091634267285,
                0.095091634267285,
                0.095091634267285,
                0.103217370534718,
                0.103217370534718,
                0.103217370534718,
                0.032458497623198,
                0.032458497623198,
                0.032458497623198,
                0.027230314174435,
                0.027230314174435,
                0.027230314174435,
                0.027230314174435,
                0.027230314174435,
                0.027230314174435
            }
        }
    );

    // 19 barycentric2D quadrature, exact for polynomials up to 9 order
    rules_.insert
    (
        9,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(0.333333333333333, 0.333333333333333, 0.333333333333333),
                barycentric2D(0.020634961602525, 0.489682519198738, 0.489682519198738),
                barycentric2D(0.489682519198738, 0.020634961602525, 0.489682519198738),
                barycentric2D(0.489682519198738, 0.489682519198738, 0.020634961602525),
                barycentric2D(0.125820817014127, 0.437089591492937, 0.437089591492937),
                barycentric2D(0.437089591492937, 0.125820817014127, 0.437089591492937),
                barycentric2D(0.437089591492937, 0.437089591492937, 0.125820817014127),
                barycentric2D(0.623592928761935, 0.188203535619033, 0.188203535619033),
                barycentric2D(0.188203535619033, 0.623592928761935, 0.188203535619033),
                barycentric2D(0.188203535619033, 0.188203535619033, 0.623592928761935),
                barycentric2D(0.910540973211095, 0.044729513394453, 0.044729513394453),
                barycentric2D(0.044729513394453, 0.910540973211095, 0.044729513394453),
                barycentric2D(0.044729513394453, 0.044729513394453, 0.910540973211095),
                barycentric2D(0.036838412054736, 0.221962989160766, 0.741198598784498),
                barycentric2D(0.036838412054736, 0.741198598784498, 0.221962989160766),
                barycentric2D(0.221962989160766, 0.036838412054736, 0.741198598784498),
                barycentric2D(0.221962989160766, 0.741198598784498, 0.036838412054736),
                barycentric2D(0.741198598784498, 0.036838412054736, 0.221962989160766),
                barycentric2D(0.741198598784498, 0.221962989160766, 0.036838412054736)
            },
            List<scalar>
            {
                0.097135796282799,
                0.031334700227139,
                0.031334700227139,
                0.031334700227139,
                0.077827541004774,
                0.077827541004774,
                0.077827541004774,
                0.079647738927210,
                0.079647738927210,
                0.079647738927210,
                0.025577675658698,
                0.025577675658698,
                0.025577675658698,
                0.043283539377289,
                0.043283539377289,
                0.043283539377289,
                0.043283539377289,
                0.043283539377289,
                0.043283539377289
            }
        }
    );


    // 25 barycentric2D quadrature, exact for polynomials up to 10 order
    rules_.insert
    (
        10,
        quadratureRule
        {
            List<barycentric2D>
            {
                barycentric2D(0.333333333333333, 0.333333333333333, 0.333333333333333),
                barycentric2D(0.028844733232685, 0.485577633383657, 0.485577633383657),
                barycentric2D(0.485577633383657, 0.028844733232685, 0.485577633383657),
                barycentric2D(0.485577633383657, 0.485577633383657, 0.028844733232685),
                barycentric2D(0.781036849029926, 0.109481575485037, 0.109481575485037),
                barycentric2D(0.109481575485037, 0.781036849029926, 0.109481575485037),
                barycentric2D(0.109481575485037, 0.109481575485037, 0.781036849029926),
                barycentric2D(0.141707219414880, 0.307939838764121, 0.550352941820999),
                barycentric2D(0.141707219414880, 0.550352941820999, 0.307939838764121),
                barycentric2D(0.307939838764121, 0.141707219414880, 0.550352941820999),
                barycentric2D(0.307939838764121, 0.550352941820999, 0.141707219414880),
                barycentric2D(0.550352941820999, 0.141707219414880, 0.307939838764121),
                barycentric2D(0.550352941820999, 0.307939838764121, 0.141707219414880),
                barycentric2D(0.025003534762686, 0.246672560639903, 0.728323904597411),
                barycentric2D(0.025003534762686, 0.728323904597411, 0.246672560639903),
                barycentric2D(0.246672560639903, 0.025003534762686, 0.728323904597411),
                barycentric2D(0.246672560639903, 0.728323904597411, 0.025003534762686),
                barycentric2D(0.728323904597411, 0.025003534762686, 0.246672560639903),
                barycentric2D(0.728323904597411, 0.246672560639903, 0.025003534762686),
                barycentric2D(0.009540815400299, 0.066803251012200, 0.923655933587500),
                barycentric2D(0.009540815400299, 0.923655933587500, 0.066803251012200),
                barycentric2D(0.066803251012200, 0.009540815400299, 0.923655933587500),
                barycentric2D(0.066803251012200, 0.923655933587500, 0.009540815400299),
                barycentric2D(0.923655933587500, 0.009540815400299, 0.066803251012200),
                barycentric2D(0.923655933587500, 0.066803251012200, 0.009540815400299)
            },
            List<scalar>
            {
                0.090817990382754,
                0.036725957756467,
                0.036725957756467,
                0.036725957756467,
                0.045321059435528,
                0.045321059435528,
                0.045321059435528,
                0.072757916845420,
                0.072757916845420,
                0.072757916845420,
                0.072757916845420,
                0.072757916845420,
                0.072757916845420,
                0.028327242531057,
                0.028327242531057,
                0.028327242531057,
                0.028327242531057,
                0.028327242531057,
                0.028327242531057,
                0.009421666963733,
                0.009421666963733,
                0.009421666963733,
                0.009421666963733,
                0.009421666963733,
                0.009421666963733
            }
        }
    );}



const Map<triQuadrature::quadratureRule>& triQuadrature::rules()
{
    if (rules_.size() == 0)
    {
        constructRules();
    }

    return rules_;
}


label triQuadrature::chooseOrder(const label requestedOrder)
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

tmp<Field<point>> triQuadrature::barycentricToPoint
(
    const List<barycentric2D>& localPts
) const
{
    tmp<Field<point>> tglobalPts(new Field<point>(localPts.size()));
    Field<point>& globalPts = tglobalPts.ref();

    forAll(globalPts, pointI)
    {
        globalPts[pointI] =
            localPts[pointI].a()*tri_.a()
          + localPts[pointI].b()*tri_.b()
          + localPts[pointI].c()*tri_.c();
    }

    return tglobalPts;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


triQuadrature::triQuadrature
(
    const point& a,
    const point& b,
    const point& c,
    const label& order
)
:
    quadrature(order),
    tri_(a, b, c)
{
    const label chosen = chooseOrder(order);
    const auto& rule = rules()[chosen];

    weights_ = rule.weights;
    points_ = barycentricToPoint(rule.points);
}


triQuadrature::triQuadrature
(
    const triangle<point, const point&>& tri,
    const label& order
)
:
    quadrature(order),
    tri_(tri.a(), tri.b(), tri.c())
{
    const label chosen = chooseOrder(order);
    const auto& rule = rules()[chosen];

    weights_ = rule.weights;
    points_ = barycentricToPoint(rule.points);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const List<point>& triQuadrature::points() const
{
    return points_;
}


const List<scalar>& triQuadrature::weights() const
{
    return weights_;
}


label triQuadrature::nPoints() const
{
    return points_.size();
}


label triQuadrature::nPoints(label order)
{
    const auto& rule = rules()[chooseOrder(order)];
    return rule.points.size();
}

} // End namespace Foam

// ************************************************************************* //
