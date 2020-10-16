/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "data.H"
#include "Time.H"
#include "solverPerformance.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::data::setSolverPerformance
(
    const word& name,
    const SolverPerformance<Type>& sp
) const
{
    dictionary& dict = const_cast<dictionary&>(solverPerformanceDict());

    List<SolverPerformance<Type>> perfs;

    if (prevTimeIndex_ != this->time().timeIndex())
    {
        // Reset solver performance between iterations
        prevTimeIndex_ = this->time().timeIndex();
        dict.clear();
    }
    else
    {
        dict.readIfPresent(name, perfs);
    }

    // Following what is done in foam-extend-4.0 solution.C: 
    // "Only the first iteration and the current iteration residuals are
    // required, so the current iteration residual replaces the previous one and
    // only the first iteration is always present, VS 2017-11-27"
    if (perfs.size() < 2)
    {
        // Append to list
        perfs.setSize(perfs.size() + 1, sp);
    }
    else
    {
        perfs.last() = sp;
    }

    dict.set(name, perfs);
}


template<class Type>
void Foam::data::setSolverPerformance
(
    const SolverPerformance<Type>& sp
) const
{
    setSolverPerformance(sp.fieldName(), sp);
}


// ************************************************************************* //
