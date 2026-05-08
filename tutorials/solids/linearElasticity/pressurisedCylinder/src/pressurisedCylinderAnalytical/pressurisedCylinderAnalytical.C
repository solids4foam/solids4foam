/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    solids4foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    solids4foam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "pressurisedCylinderAnalytical.H"
#include "compatibilityFunctions.H"

namespace Foam
{

tensor cylindricalToCartesianRotation(const point& pt)
{
    const Foam::scalar r = Foam::mag(Foam::vector(pt.x(), pt.y(), 0));

    if (r < Foam::SMALL)
    {
        return Foam::tensor
        (
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        );
    }

    const Foam::scalar c = pt.x()/r;
    const Foam::scalar s = pt.y()/r;

    return Foam::tensor
    (
        c, -s, 0,
        s,  c, 0,
        0,  0, 1
    );
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


 defineTypeNameAndDebug(pressurisedCylinderAnalytical, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


pressurisedCylinderAnalytical::pressurisedCylinderAnalytical()
:
    Ri_(0),
    Ro_(0),
    p_(0),
    E_(0),
    nu_(0),
    A_(0)
{}


pressurisedCylinderAnalytical::pressurisedCylinderAnalytical
(
    const pressurisedCylinderAnalytical& pc
)
:
    Ri_(pc.Ri_),
    Ro_(pc.Ro_),
    p_(pc.p_),
    E_(pc.E_),
    nu_(pc.nu_),
    A_(p_*sqr(Ri_)/(sqr(Ro_) - sqr(Ri_)))
{}


pressurisedCylinderAnalytical::pressurisedCylinderAnalytical
(
    const dictionary& dict
)
:
    Ri_(readScalar(dict.lookup("Ri"))),
    Ro_(readScalar(dict.lookup("Ro"))),
    p_(readScalar(dict.lookup("p"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    A_(p_*sqr(Ri_)/(sqr(Ro_) - sqr(Ri_)))
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //


vector pressurisedCylinderAnalytical::displacement(const point& pt) const
{
    const vector cylindricalD = cylindricalDisplacement(pt);
    const tensor R = cylindricalToCartesianRotation(pt);

    return R & cylindricalD;
}


vector pressurisedCylinderAnalytical::cylindricalDisplacement
(
    const point& pt
) const
{
    if (E_ < SMALL || Ro_ < SMALL || Ri_ < SMALL || p_ < SMALL)
    {
        FatalErrorInFunction
            << "Ro, Ri, p and E must be greater than 0!" << exit(FatalError);
    }

    if (Ro_ < Ri_)
    {
        FatalErrorInFunction
            << "Ro should be greater than Ri!" << exit(FatalError);
    }

    const scalar r = mag(vector(pt.x(), pt.y(), 0));

    const scalar ur =
        (p_/E_)*sqr(Ri_)/(sqr(Ro_) - sqr(Ri_))
       *((1.0 - nu_)*r + (1.0 + nu_)*sqr(Ro_)/r);

    if (r < SMALL)
    {
        return vector::zero;
    }

    return vector(ur, 0.0, 0.0);
}


tmp<vectorField> pressurisedCylinderAnalytical::displacement
(
    const vectorField& locations
) const
{
    // Prepare the result field
    tmp<vectorField> tresult(new vectorField(locations.size(), vector::zero));
    vectorField& result = tmpRef(tresult);

    // Calculate the result per location
    forAll(result, i)
    {
        result[i] = displacement(locations[i]);
    }

    return tresult;
}


symmTensor pressurisedCylinderAnalytical::stress(const point& pt) const
{
    const symmTensor cylindricalS = cylindricalStress(pt);
    const tensor R = cylindricalToCartesianRotation(pt);

    return symm(R & (cylindricalS & R.T()));
}


symmTensor pressurisedCylinderAnalytical::cylindricalStress
(
    const point& pt
) const
{
    if (E_ < SMALL || Ro_ < SMALL || Ri_ < SMALL || p_ < SMALL)
    {
        FatalErrorInFunction
            << "Ro, Ri, p and E must be greater than 0!" << exit(FatalError);
    }

    if (Ro_ < Ri_)
    {
        FatalErrorInFunction
            << "Ro should be greater than Ri!" << exit(FatalError);
    }

    const scalar r = mag(vector(pt.x(), pt.y(), 0));

    if (r < SMALL)
    {
        return symmTensor::zero;
    }

    const scalar sigR = A_*(1.0 - sqr(Ro_)/sqr(r));
    const scalar sigT = A_*(1.0 + sqr(Ro_)/sqr(r));

    return symmTensor
    (
        sigR,
        0.0,
        0.0,
        sigT,
        0.0,
        0.0
    );
}


void Foam::pressurisedCylinderAnalytical::write(Ostream& os) const
{
    os.writeKeyword("Ri")
        << Ri_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ro")
        << Ro_ << token::END_STATEMENT << nl;
    os.writeKeyword("E")
        << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("p")
        << p_ << token::END_STATEMENT << nl;
    os.writeKeyword("nu")
        << nu_ << token::END_STATEMENT << nl;
}

}  // End namespace Foam

// ************************************************************************* //
