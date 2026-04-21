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
    setPlateHoleBC

Description
    FunctionObject updates displacement boundary conditions for stress analysis

Author
    Zeljko Tukovic, FSB.  All rights reserved

SourceFiles
    setPlateHoleBC.C

\*---------------------------------------------------------------------------*/

#ifndef setPlateHoleBC_H
#define setPlateHoleBC_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

#include "cylindricalCS.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setPlateHoleBC Declaration
\*---------------------------------------------------------------------------*/

class setPlateHoleBC
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Hole radius
        scalar a_;

        //- Loading traction
        scalar T_;

    // Private Member Functions

        //- Check if this is plane stress simulation
        bool planeStress() const;

        //- Plate hole solution - stress
        symmTensor plateHoleStress(const vector& C) const;

        //- Plate hole solution - stress
        scalar plateHoleHydPressure(const vector& C) const;

        //- Plate hole solution - displacement
        vector plateHoleDisplacement(const vector& C) const;

        //- Disallow default bitwise copy construct
        setPlateHoleBC
        (
            const setPlateHoleBC&
        );

        //- Disallow default bitwise assignment
        void operator=(const setPlateHoleBC&);

        //- Update boundary condition
        void setBC();

        //- Calculate error
        void calcError() const;

public:

    //- Runtime type information
    TypeName("setPlateHoleBC");


    // Constructors

        //- Construct from components
        setPlateHoleBC
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
// #if FOAMEXTEND > 40
        virtual bool execute(const bool forceWrite);
// #else
//         virtual bool execute();
// #endif

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
