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
    Foam::mmsFsiErrorsFunctionObject

Description
    Error norms of the manufactured FSI solution, written every time step to
    postProcessing/mmsFsiErrors/errors.dat. The L1, L2 and LInf norms of the
    error magnitude are given for:

    - U, p: the fluid velocity and kinematic pressure at the cell centres of
      the current (moving) fluid mesh, volume weighted; pMean_L2 is the L2
      norm of the pressure error after removing its mean;
    - D, Us: the solid displacement and velocity at the cell centres of the
      undeformed solid mesh, volume weighted; the velocity is the ddt of D
      with the case's ddt scheme;
    - intD: the displacement of the solid interface faces, area weighted;
    - intTrac: the traction applied to the solid interface by the
      solidTraction condition (traction - pressure*n), against the exact
      Cauchy traction, area weighted;
    - intP: the fluid pressure on the fluid interface faces, area weighted.

    With writeFields true, the fluid errors UError and pError are also
    written at the write times.

Usage
    \verbatim
    mmsErrors
    {
        type            mmsFsiErrors;
        libs            (mmsFsi);

        // Optional entries (defaults shown)
        fluidRegion     fluid;
        solidRegion     solid;
        interfacePatch  interface;
        writeFields     false;
    }
    \endverbatim

SourceFiles
    mmsFsiErrorsFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef mmsFsiErrorsFunctionObject_H
#define mmsFsiErrorsFunctionObject_H

#include "functionObject.H"
#include "OFstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

/*---------------------------------------------------------------------------*\
                 Class mmsFsiErrorsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class mmsFsiErrorsFunctionObject
:
    public functionObject
{
    // Private Data

        //- Reference to the main Time
        const Time& time_;

        //- Name of the fluid region
        word fluidRegion_;

        //- Name of the solid region
        word solidRegion_;

        //- Name of the interface patch in both regions
        word interfacePatch_;

        //- Write the fluid error fields at the write times
        bool writeFields_;

        //- Output file
        autoPtr<OFstream> filePtr_;


    // Private Member Functions

        //- Calculate and write the errors
        bool writeErrors();


public:

    //- Runtime type information
    TypeName("mmsFsiErrors");


    // Constructors

        //- Construct from Time and dictionary
        mmsFsiErrorsFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        mmsFsiErrorsFunctionObject(const mmsFsiErrorsFunctionObject&) =
            delete;

        //- No copy assignment
        void operator=(const mmsFsiErrorsFunctionObject&) = delete;


    //- Destructor
    virtual ~mmsFsiErrorsFunctionObject() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Calculate and write the errors
        virtual bool execute();

        //- Do nothing: the errors are written by execute
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
