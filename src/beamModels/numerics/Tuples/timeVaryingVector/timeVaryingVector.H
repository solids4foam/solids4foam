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

Class
    Foam::timeVaryingVector

Description

SourceFiles
    timeVaryingVector.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingVector_H
#define timeVaryingVector_H

#include "vector.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
// class Ostream;

// Forward declaration of friend functions and operators
class timeVaryingVector;
Istream& operator>>(Istream&, timeVaryingVector&);
// Ostream& operator<<(Ostream&, const timeVaryingVector&);


/*---------------------------------------------------------------------------*\
                          Class timeVaryingVector Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingVector
:
    public vector
{
    // Private data

        //- Vector time series
        interpolationTable<vector> timeSeries_;


    // Private Member Functions

        // //- Disallow default bitwise copy construct
        // timeVaryingVector(const timeVaryingVector&);

        // //- Disallow default bitwise assignment
        // void operator=(const timeVaryingVector&);


public:

    // Static data members

        // //- Static data staticData
        // static const dataType staticData;


    // Constructors

        //- Construct null
        timeVaryingVector();

        //- Construct from components
        timeVaryingVector(const vector& v);

        //- Construct from Istream
        timeVaryingVector(Istream&);

        //- Construct given three components
        timeVaryingVector(const scalar& vx, const scalar& vy, const scalar& vz);


    // Selectors

        // //- Select null constructed
        // static autoPtr<timeVaryingVector> New();


    //- Destructor
    ~timeVaryingVector();


    // Member Functions

        // Access

        // Returne reference to interpolation table
        interpolationTable<vector>& timeSeries()
        {
            return timeSeries_;
        }

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const timeVaryingVector&);

        vector operator()(const scalar t) const
        {
            return timeSeries_(t);
        }


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        // friend Istream& operator>>(Istream&, timeVaryingVector&);
        // friend Ostream& operator<<(Ostream&, const timeVaryingVector&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "timeVaryingVectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
