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

#include "cantileverAnalytical.H"
#include "compatibilityFunctions.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cantileverAnalytical, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cantileverAnalytical::cantileverAnalytical()
:
    P_(0),
    E_(0),
    nu_(0),
    L_(0),
    D_(0),
    I_(0)
{}


Foam::cantileverAnalytical::cantileverAnalytical(const cantileverAnalytical& ca)
:
    P_(ca.P_),
    E_(ca.E_),
    nu_(ca.nu_),
    L_(ca.L_),
    D_(ca.D_),
    I_(ca.I_)
{}


Foam::cantileverAnalytical::cantileverAnalytical(const dictionary& dict)
:
    P_(readScalar(dict.lookup("P"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    L_(readScalar(dict.lookup("L"))),
    D_(readScalar(dict.lookup("D"))),
    I_(Foam::pow(D_, 3.0)/12.0)
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //


Foam::vector Foam::cantileverAnalytical::displacement
(
    const vector& location
) const
{
    if (E_ < SMALL || L_ < SMALL || D_ < SMALL || I_ < SMALL)
    {
        FatalErrorInFunction
            << "E, L, D and I must be greater than 0!" << exit(FatalError);
    }

    // Extract x and y location components
    const scalar x(location.component(vector::X));
    const scalar y(location.component(vector::Y));

    // Return the displacement per Augarde and Deeks
    return
        vector
        (
            (P_*y/(6*E_*I_))*((6*L_ - 3*x)*x + (2 + nu_)*(y*y - D_*D_/4.0)),
           -(P_/(6*E_*I_))
           *(
                3*nu_*y*y*(L_ - x) + (4 + 5*nu_)*D_*D_*x/4 + (3*L_ - x)*x*x
            ),
            0
        );
}


Foam::tmp<Foam::vectorField> Foam::cantileverAnalytical::displacement
(
    const vectorField& locations
) const
{
    if (E_ < SMALL || L_ < SMALL || D_ < SMALL || I_ < SMALL)
    {
        FatalErrorInFunction
            << "E, L, D and I must be greater than 0!" << exit(FatalError);
    }

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


Foam::vector Foam::cantileverAnalytical::traction
(
    const vector& location
) const
{
    if (E_ < SMALL || L_ < SMALL || D_ < SMALL || I_ < SMALL)
    {
        FatalErrorInFunction
            << "E, L, D and I must be greater than 0!" << exit(FatalError);
    }

    // Extract x and y location components
    const scalar x(location.component(vector::X));
    const scalar y(location.component(vector::Y));

    // Return the traction on the x plane
    return
        vector
        (
            P_*(L_ - x)*y/I_,
            -(P_/(2*I_))*(D_*D_/4 - y*y),
            0
        );
}


Foam::tmp<Foam::vectorField> Foam::cantileverAnalytical::traction
(
    const vectorField& locations
) const
{
    if (E_ < SMALL || L_ < SMALL || D_ < SMALL || I_ < SMALL)
    {
        FatalErrorInFunction
            << "E, L, D and I must be greater than 0!" << exit(FatalError);
    }

    // Prepare the result field
    tmp<vectorField> tresult(new vectorField(locations.size(), vector::zero));
    vectorField& result = tmpRef(tresult);

    // Calculate the result per location
    forAll(result, i)
    {
        result[i] = traction(locations[i]);
    }

    return tresult;
}


Foam::symmTensor Foam::cantileverAnalytical::stress
(
    const vector& location
) const
{
    if (E_ < SMALL || L_ < SMALL || D_ < SMALL || I_ < SMALL)
    {
        FatalErrorInFunction
            << "E, L, D and I must be greater than 0!" << exit(FatalError);
    }

    // Calculate the traction on the x plane
    const vector trac(traction(location));

    // Return the stress tensor
    return
        symmTensor
        (
            trac.x(), trac.y(), 0.0,
                      0.0,      0.0,
                                0.0
        );
}


void Foam::cantileverAnalytical::write(Ostream& os) const
{
    os.writeKeyword("P")
        << P_ << token::END_STATEMENT << nl;
    os.writeKeyword("E")
        << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("nu")
        << nu_ << token::END_STATEMENT << nl;
    os.writeKeyword("L")
        << L_ << token::END_STATEMENT << nl;
    os.writeKeyword("D")
        << D_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
