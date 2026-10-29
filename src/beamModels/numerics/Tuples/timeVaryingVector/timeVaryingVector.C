/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "timeVaryingVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// const dataType Foam::timeVaryingVector::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingVector::timeVaryingVector()
:
    vector(),
    timeSeries_()
{
    // Info << "Creating zero vector" << endl;
}


Foam::timeVaryingVector::timeVaryingVector(const vector& v)
:
    vector(v),
    timeSeries_()
{}


Foam::timeVaryingVector::timeVaryingVector(Istream& s)
:
    vector(s),
    timeSeries_()
{
    // Info << "reading vector" << endl;
}


Foam::timeVaryingVector::timeVaryingVector
(
    const scalar& vx,
    const scalar& vy,
    const scalar& vz
)
:
    vector(vx, vy, vz),
    timeSeries_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Foam::autoPtr<Foam::timeVaryingVector> Foam::timeVaryingVector::New()
// {
//     return autoPtr<timeVaryingVector>(new timeVaryingVector);
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeVaryingVector::~timeVaryingVector()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::timeVaryingVector::operator=(const timeVaryingVector& tvv)
{
    Info << "using = operator" << endl;
}



// void Foam::timeVaryingVector::operator=(const timeVaryingVector& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("Foam::timeVaryingVector::operator=(const Foam::timeVaryingVector&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// }

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


Foam::Istream& Foam::operator>>
(
    Istream& is,
    timeVaryingVector& tvv
)
{
    fileName forceFileName;
    //bounds::repeatableBounding::CLAMP;

    is >> forceFileName;

    tvv.timeSeries() = interpolationTable<vector>(forceFileName); //AT : bounds::repeatableBounds::CLAMP, needs to be implemented.
    //tvv.timeSeries().outOfBounds(interpolationTable<vector>::CLAMP);

    // Info << forceFileName << endl;

    vector& v = tvv;
    v = tvv.timeSeries()(0);

    // is >> v;

    return is;
}


// ************************************************************************* //
