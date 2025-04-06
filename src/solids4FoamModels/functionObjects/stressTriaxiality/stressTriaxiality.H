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
    stressTriaxiality

Description
    FunctionObject that outputs the stress triaxiality field:
    where: stressTriaxiality = hydrostaticStress/equivalentStress;

    @verbatim
    stressTriaxiality1
    {
        // Type of functionObject
        type stressTriaxiality;

        // Where to load it from (if not already in solver)
        // Note: this is not required if the solver already loads this library.
        libs ("libsolids4FoamModels.so");

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;
    );
    @endverbatim

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    stressTriaxiality.C

\*---------------------------------------------------------------------------*/

#ifndef stressTriaxiality_H
#define stressTriaxiality_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class stressTriaxiality Declaration
\*---------------------------------------------------------------------------*/

class stressTriaxiality
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& runTime_;

        //- Reference to mesh region
        const fvMesh& mesh_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        stressTriaxiality
        (
            const stressTriaxiality&
        );

        //- Disallow default bitwise assignment
        void operator=(const stressTriaxiality&);


public:

    //- Runtime type information
    TypeName("stressTriaxiality");


    // Constructors

        //- Construct from components
        stressTriaxiality
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND
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
