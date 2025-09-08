/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
    Copyright (C) 2019 Simone Bna
    Copyright (C) 2020 Stefano Zampini
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

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"
#include "error.H"
#include "petscUtils.H"
#include "petscErrorHandling.H"

// For older PETSc. Not strictly correct for 64-bit compilations,
// but adequate for transitional code
#ifndef PetscInt_FMT
#define PetscInt_FMT "D"
#endif

// For PETSc < 3.19.0
#ifndef PETSC_SUCCESS
#define PETSC_SUCCESS 0
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::PetscUtils::setFlags
(
    const word& prefix,
    const dictionary& dict,
    const bool verbose
)
{
    for (const entry& e : dict)
    {
        if (e.isDict())
        {
            // Skip sub-dicts
            WarningInFunction
                << "Skipping sub-dict " << e.keyword() << " in " << dict.name()
                << endl;

            continue;
        }

        const word key = '-' + prefix + e.keyword();

        const ITstream& is = e.stream();
        const label nTok  = is.size();

        if (nTok == 0)
        {
            // No value → PETSc flag
            if (verbose)
            {
                Info<< key << nl;
            }
            AssertPETSc(PetscOptionsSetValue(NULL, key.c_str(), nullptr));
        }
        else
        {
            // Has a value (assume single-token PETSc option values, e.g. "newtonls")
            const word val = e.get<word>();  // safe because nTok>0
            if (verbose)
            {
                Info<< key << ' ' << val << nl;
            }
            AssertPETSc(PetscOptionsSetValue(NULL, key.c_str(), val.c_str()));
        }
    }
}


// ************************************************************************* //
