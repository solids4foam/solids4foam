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
    squarePlateAnalyticalSolution

Description
    FunctionObject that generates the analytical deflection for the thin square
    plate with uniform transverse pressure.

    The solution is taken from:
    Timoshenko, S., & Woinowsky-Krieger, S. (1959). Theory of plates and shells.

    It is assumed that the origin is at the centre of the plate.
    The plate dimensions are a x b, thickness of the plate is h.

Author
    Ivan Batistic & Philip Cardiff, UCD. All rights reserved.

SourceFiles
    squarePlateAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef squarePlateAnalyticalSolution_H
#define squarePlateAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class squarePlateAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class squarePlateAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Applied transverse pressure
        const scalar p_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Length of the plate (in x direction)
        const scalar a_;

        //- Width of the plate (in y direction)
        const scalar b_;

        //- Thickness of plate
        const scalar h_;

        //- Boundary conditions combination
        const word boundaryConditions_;

        //- Calculate cell displacement fields
        const Switch cellDisplacement_;

        //- Calculate point displacement fields
        const Switch pointDisplacement_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        squarePlateAnalyticalSolution
        (
            const squarePlateAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const squarePlateAnalyticalSolution&);


public:

    //- Runtime type information
    TypeName("squarePlateAnalyticalSolution");


    // Constructors

        //- Construct from components
        squarePlateAnalyticalSolution
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#ifdef FOAMEXTEND
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

#ifdef OPENFOAM_NOT_EXTEND
        //- Write
        virtual bool write();
#endif

#ifndef OPENFOAM_NOT_EXTEND
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
