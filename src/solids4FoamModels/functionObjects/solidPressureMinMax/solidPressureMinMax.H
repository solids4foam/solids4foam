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
    solidPressureMinMax

Description
    FunctionObject that writes the global min and max of the hydrostatic
    pressure in the solid at every time step.

    The hydrostatic pressure is defined as

        p = -tr(sigma) / 3

    where sigma is the Cauchy stress tensor.  Both min and max are written to
    a plain text data file:

        postProcessing/<startTime>/solidPressureMinMax.dat

    File columns:
        Time  p_min  p_max  trSigma_min  trSigma_max

    Dictionary example:

        pressureMinMax
        {
            type    solidPressureMinMax;
        }

Author
    Philip Cardiff, UCD.

SourceFiles
    solidPressureMinMax.C

\*---------------------------------------------------------------------------*/

#ifndef solidPressureMinMax_H
#define solidPressureMinMax_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidPressureMinMax Declaration
\*---------------------------------------------------------------------------*/

class solidPressureMinMax
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- History file ptr
        autoPtr<OFstream> historyFilePtr_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        solidPressureMinMax(const solidPressureMinMax&);

        //- Disallow default bitwise assignment
        void operator=(const solidPressureMinMax&);


public:

    //- Runtime type information
    TypeName("solidPressureMinMax");


    // Constructors

        //- Construct from components
        solidPressureMinMax
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
