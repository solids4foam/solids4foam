/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mechanicalConstitutiveLawResponse

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
//  This enum specifies which (if any) tangent quantities a mechanical
//  constitutive law should compute during evaluation. The request is
//  advisory: a law may choose to ignore unsupported requests.
enum class tangentRequest
{
    none,       //- Stress only
    scalar      //- Scalar approximate tangent
    // future:
    // symmTensor,
    // fourthOrder
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
