/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solution.H"
#include "objectRegistry.H"

// These are for old syntax compatibility:
#include "BICCG.H"
#include "ICCG.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::debug::debugSwitch
Foam::solution::debug
(
    "solution",
    0
);

// List of sub-dictionaries to rewrite
//! @cond localScope
static const Foam::List<Foam::word> subDictNames
(
    Foam::IStringStream("(preconditioner smoother)")()
);
//! @endcond localScope


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solution::read(const dictionary& dict)
{
    if (dict.found("cache"))
    {
        cache_ = dict.subDict("cache");
        caching_ = cache_.lookupOrDefault("active", true);
    }

    if (dict.found("relaxationFactors"))
    {
        const dictionary& relaxDict(dict.subDict("relaxationFactors"));
        if (relaxDict.found("fields") || relaxDict.found("equations"))
        {
            if (relaxDict.found("fields"))
            {
                fieldRelaxDict_ = relaxDict.subDict("fields");
            }

            if (relaxDict.found("equations"))
            {
                eqnRelaxDict_ = relaxDict.subDict("equations");
            }
        }
        else
        {
            // backwards compatibility
            fieldRelaxDict_.clear();

            const wordList entryNames(relaxDict.toc());
            forAll(entryNames, i)
            {
                const word& e = entryNames[i];
                scalar value = readScalar(relaxDict.lookup(e));

                if (e(0, 1) == "p")
                {
                    fieldRelaxDict_.add(e, value);
                }
                else if (e.length() >= 3)
                {
                    if (e(0, 3) == "rho")
                    {
                        fieldRelaxDict_.add(e, value);
                    }
                }

            }

            eqnRelaxDict_ = relaxDict;
        }

        fieldRelaxDefault_ =
            fieldRelaxDict_.lookupOrDefault<scalar>("default", 0.0);

        eqnRelaxDefault_ =
            eqnRelaxDict_.lookupOrDefault<scalar>("default", 0.0);

        if (debug)
        {
            Info<< "relaxation factors:" << nl
                << "fields: " << fieldRelaxDict_ << nl
                << "equations: " << eqnRelaxDict_ << endl;
        }
    }


    if (dict.found("solvers"))
    {
        solvers_ = dict.subDict("solvers");
        upgradeSolverDict(solvers_);
    }

    if (dict.found("solverPerformance"))
    {
        solverPerformance_ = dict.subDict("solverPerformance");
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution(const objectRegistry& obr, const fileName& dictName)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            IOobject::READ_IF_PRESENT,  // Allow default dictionary creation
            IOobject::NO_WRITE
        )
    ),
    cache_(dictionary::null),
    caching_(false),
    fieldRelaxDict_(dictionary::null),
    eqnRelaxDict_(dictionary::null),
    fieldRelaxDefault_(0),
    eqnRelaxDefault_(0),
    solvers_(dictionary::null),
    solverPerformance_(dictionary::null),
    prevTimeIndex_(0)
{
    if (!headerOk())
    {
        if (debug)
        {
            InfoIn
            (
                "Foam::solution::solution(const objectRegistry& obr, "
                "const fileName& dictName)"
            )   << "Solution dictionary not found.  Adding default entries"
                << endl;
        }
    }

    read(solutionDict());

    set("solverPerformance", dictionary());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::solution::upgradeSolverDict
(
    dictionary& dict,
    const bool verbose
)
{
    label nChanged = 0;

    // backward compatibility:
    // recast primitive entries into dictionary entries
    forAllIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            Istream& is = iter().stream();
            word name(is);
            dictionary subdict;

            if (name == "BICCG")
            {
                // special treatment for very old syntax
                subdict = BICCG::solverDict(is);
            }
            else if (name == "ICCG")
            {
                // special treatment for very old syntax
                subdict = ICCG::solverDict(is);
            }
            else
            {
                subdict.add("solver", name);
                subdict <<= dictionary(is);

                // preconditioner and smoother entries can be
                // 1) primitiveEntry w/o settings,
                // 2) or a dictionaryEntry.
                // transform primitiveEntry with settings -> dictionaryEntry
                forAll(subDictNames, dictI)
                {
                    const word& dictName = subDictNames[dictI];
                    entry* ePtr = subdict.lookupEntryPtr(dictName,false,false);

                    if (ePtr && !ePtr->isDict())
                    {
                        Istream& is = ePtr->stream();
                        is >> name;

                        if (!is.eof())
                        {
                            dictionary newDict;
                            newDict.add(dictName, name);
                            newDict <<= dictionary(is);

                            subdict.set(dictName, newDict);
                        }
                    }
                }
            }


            // write out information to help people adjust to the new syntax
            if (verbose && Pstream::master())
            {
                Info<< "// using new solver syntax:\n"
                    << iter().keyword() << subdict << endl;
            }

            // overwrite with dictionary entry
            dict.set(iter().keyword(), subdict);

            nChanged++;
        }
    }

    return nChanged;
}


bool Foam::solution::cache(const word& name) const
{
    if (caching_)
    {
        if (debug)
        {
            Info<< "Cache: find entry for " << name
                << ": " << Switch(cache_.found(name)) << endl;
        }

        return cache_.found(name);
    }
    else
    {
        return false;
    }
}


bool Foam::solution::relaxField(const word& name) const
{
    if (debug)
    {
        Info<< "Field relaxation factor for " << name
            << " is " << (fieldRelaxDict_.found(name) ? "set" : "unset")
            << endl;
    }

    return fieldRelaxDict_.found(name) || fieldRelaxDict_.found("default");
}


bool Foam::solution::relaxEquation(const word& name) const
{
    if (debug)
    {
        Info<< "Find equation relaxation factor for " << name << endl;
    }

    return eqnRelaxDict_.found(name) || eqnRelaxDict_.found("default");
}


Foam::scalar Foam::solution::fieldRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup variable relaxation factor for " << name << endl;
    }

    if (fieldRelaxDict_.found(name))
    {
        return readScalar(fieldRelaxDict_.lookup(name));
    }
    else if (fieldRelaxDefault_ > SMALL)
    {
        return fieldRelaxDefault_;
    }
    else
    {
        FatalIOErrorIn
        (
            "Foam::solution::fieldRelaxationFactor(const word&)",
            fieldRelaxDict_
        )   << "Cannot find variable relaxation factor for '" << name
            << "' or a suitable default value."
            << exit(FatalIOError);

        return 0;
    }
}


Foam::scalar Foam::solution::equationRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup equation relaxation factor for " << name << endl;
    }

    if (eqnRelaxDict_.found(name))
    {
        return readScalar(eqnRelaxDict_.lookup(name));
    }
    else if (eqnRelaxDefault_ > SMALL)
    {
        return eqnRelaxDefault_;
    }
    else
    {
        FatalIOErrorIn
        (
            "Foam::solution::eqnRelaxationFactor(const word&)",
            eqnRelaxDict_
        )   << "Cannot find equation relaxation factor for '" << name
            << "' or a suitable default value."
            << exit(FatalIOError);

        return 0;
    }
}


const Foam::dictionary& Foam::solution::solutionDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solverDict(const word&)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


const Foam::dictionary& Foam::solution::solver(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solver(const word&)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        read(solutionDict());

        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::solution::writeData(Ostream& os) const
{
    // Write direct entries of the solution dictionary
    // HJ, 16/Feb/2010
    dictionary::write(os, false);

    return true;
}

Foam::dictionary& Foam::solution::solverPerformanceDict() const
{
    return solverPerformance_;
}


void Foam::solution::setSolverPerformance
(
    const word& name,
    const lduSolverPerformance& sp
) const
{
    List<lduSolverPerformance> perfs;

    if (prevTimeIndex_ != this->time().timeIndex())
    {
        // Reset solver performance between iterations
        prevTimeIndex_ = this->time().timeIndex();
        solverPerformance_.clear();
    }
    else
    {
        solverPerformance_.readIfPresent(name, perfs);
    }

    // Only the first iteration and the current iteration residuals are
    // required, so the current iteration residual replaces the previous one and
    // only the first iteration is always present, VS 2017-11-27
    if (perfs.size() < 2)
    {
        // Append to list
        perfs.setSize(perfs.size() + 1, sp);
    }
    else
    {
        perfs.last() = sp;
    }

    solverPerformance_.set(name, perfs);
}


void Foam::solution::setSolverPerformance
(
    const lduSolverPerformance& sp
) const
{
    setSolverPerformance(sp.fieldName(), sp);
}

// ************************************************************************* //
