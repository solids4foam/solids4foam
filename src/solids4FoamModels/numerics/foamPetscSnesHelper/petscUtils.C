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

#ifdef USE_PETSC

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
#ifdef OPENFOAM_COM
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
#else // FOAMEXTEND and OPENFOAM_ORG
    // foam-extend: iterate with FOAM's list iterator macro
    forAllConstIter(dictionary, dict, iter)
    {
        const entry& e = iter();

        if (e.isDict())
        {
            WarningInFunction
                << "Skipping sub-dict " << e.keyword() << " in " << dict.name()
                << endl;
            continue;
        }

        const word key = '-' + prefix + e.keyword();

        // In foam-extend, use the token stream
        if (!e.isStream() || e.stream().size() == 0)
        {
            // No value -> PETSc flag
            if (verbose) Info<< key << nl;
            AssertPETSc(PetscOptionsSetValue(NULL, key.c_str(), nullptr));
        }
        else
        {
            // Read the first token as the value
            // stream API is non-const
            ITstream& is = const_cast<ITstream&>(e.stream());
            word val;
            is >> val;   // assume single-token PETSc option values
            if (verbose) Info<< key << ' ' << val << nl;
            AssertPETSc(PetscOptionsSetValue(NULL, key.c_str(), val.c_str()));
        }
    }
#endif // FOAMEXTEND
}

// ************************************************************************* //

#endif // #ifdef USE_PETSC
