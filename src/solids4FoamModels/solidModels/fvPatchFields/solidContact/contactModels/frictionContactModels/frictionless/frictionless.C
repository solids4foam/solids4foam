/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "frictionless.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(frictionless, 0);
    addToRunTimeSelectionTable(frictionContactModel, frictionless, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::frictionless::frictionless
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID
)
:
    frictionContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID
    ),
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero)
{}


// Construct as a copy
Foam::frictionless::frictionless(const frictionless& fricLaw)
:
    frictionContactModel(fricLaw),
    slaveTraction_(fricLaw.slaveTraction_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::frictionless::~frictionless()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::frictionless::autoMap(const fvPatchFieldMapper& m)
{
    frictionContactModel::autoMap(m);

    if (debug)
    {
        InfoIn
        (
            "void frictionless::autoMap(const fvPatchFieldMapper& m)"
        )   << "autoMap" << endl;
    }

#ifdef OPENFOAMFOUNDATION
    m(slaveTraction_, slaveTraction_);
#else
    slaveTraction_.autoMap(m);
#endif
}

// ************************************************************************* //
