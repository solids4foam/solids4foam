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

Class
    Foam::mechanicalConstitutiveLawTangentRequest

Description
    Enumeration defining requests for constitutive tangent information.

    TangentRequest is used to explicitly specify which (if any) constitutive
    tangent quantities should be computed during a constitutive law
    evaluation. The request expresses intent only; a mechanical constitutive
    law may choose to ignore unsupported requests.

    This mechanism allows solvers to control when tangents are evaluated
    independently of stress updates, which is particularly important for
    segregated solution strategies where stress may be recomputed many times
    while tangents are only required occasionally.

    The enumeration is intentionally minimal and extensible, allowing future
    support for higher-order or more detailed tangent representations without
    changes to the constitutive law interface.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawTangentRequest_H
#define mechanicalConstitutiveLawTangentRequest_H

namespace Foam
{

//- Requested level of constitutive tangent information.
//
//  Specifies which (if any) tangent quantities a mechanical constitutive
//  law should provide during evaluation. The request is advisory: a law
//  may ignore requests it does not support.
//  Tangent quantities are intended for use in implicit or mixed
//  formulations (e.g. Newton–Raphson schemes or mixed u–p methods).
enum class tangentRequest
{
    none,                       // Stress only; no tangent information

    scalar,                     // Scalar approximate tangent

    scalarDeviatoric,           // Deviatoric-only scalar tangent
                                // (e.g. for mixed displacement–pressure
                                // formulations)

    fourthOrder,                // Full fourth-order material tangent,
                                // represented as a 6x6 matrix in Voigt form

    fourthOrderFiniteDifference // Full fourth-order material tangent,
                                // represented as a 6x6 Voigt matrix and
                                // computed using finite differences
};

// Return the string name of a tangent request
inline const word tangentRequestName
(
    const tangentRequest r
)
{
    switch (r)
    {
        case tangentRequest::none:              return "none";
        case tangentRequest::scalar:            return "scalar";
        case tangentRequest::scalarDeviatoric:  return "scalarDeviatoric";
        case tangentRequest::fourthOrder:       return "fourthOrder";
        case tangentRequest::fourthOrderFiniteDifference:
            return "fourthOrderFiniteDifference";
    }

    return "unknown";
}


// Return the tangent request with the given name.
// Written by hand rather than with NamedEnum or Enum, whose interfaces differ
// between the supported forks
inline tangentRequest tangentRequestNamed(const word& name)
{
    if (name == "none")
    {
        return tangentRequest::none;
    }
    else if (name == "scalar")
    {
        return tangentRequest::scalar;
    }
    else if (name == "scalarDeviatoric")
    {
        return tangentRequest::scalarDeviatoric;
    }
    else if (name == "fourthOrder")
    {
        return tangentRequest::fourthOrder;
    }
    else if (name == "fourthOrderFiniteDifference")
    {
        return tangentRequest::fourthOrderFiniteDifference;
    }

    FatalErrorInFunction
        << "Unknown tangent request '" << name << "'." << nl
        << "Valid options are: none, scalar, scalarDeviatoric, fourthOrder, "
        << "fourthOrderFiniteDifference"
        << exit(FatalError);

    return tangentRequest::none;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
