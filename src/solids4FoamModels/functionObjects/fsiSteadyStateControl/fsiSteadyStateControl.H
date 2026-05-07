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
    fsiSteadyStateControl

Description
    Function object for detecting steady FSI responses from a monitored solid
    point displacement and interface force. When the monitored values change by
    less than specified tolerances over a sliding sample window, the function
    object can write the current time and end the run cleanly.

SourceFiles
    fsiSteadyStateControl.C

\*---------------------------------------------------------------------------*/

#ifndef fsiSteadyStateControl_H
#define fsiSteadyStateControl_H

#include "stateFunctionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fsiSteadyStateControl Declaration
\*---------------------------------------------------------------------------*/

class fsiSteadyStateControl
:
    public functionObjects::stateFunctionObject
{
    // Private data

        //- Name
        const word name_;

        //- Solid region name
        word solidRegionName_;

        //- Name of the forces function object to read
        word forcesObjectName_;

        //- Monitored point location
        vector point_;

        //- Displacement component to monitor
        direction displacementComponent_;

        //- Force component to monitor
        direction forceComponent_;

        //- Minimum run time before checking convergence
        scalar minTime_;

        //- Number of time samples separating convergence comparisons
        label nSamples_;

        //- Number of consecutive satisfied checks required
        label nConsecutive_;

        //- Absolute tolerance for displacement changes
        scalar displacementAbsTol_;

        //- Relative tolerance for displacement changes
        scalar displacementRelTol_;

        //- Absolute tolerance for force changes
        scalar forceAbsTol_;

        //- Relative tolerance for force changes
        scalar forceRelTol_;

        //- Lower reference magnitude for displacement relative changes
        scalar displacementRefMin_;

        //- Lower reference magnitude for force relative changes
        scalar forceRefMin_;

        //- Include viscous force contribution
        Switch includeViscousForce_;

        //- Include internal force contribution
        Switch includeInternalForce_;

        //- Write and end when steady
        Switch writeAndEnd_;

        //- Log check status
        Switch log_;

        //- Number of consecutive checks currently satisfied
        label consecutive_;

        //- Has the run already been stopped by this function object?
        bool stopped_;

        //- Displacement sample history
        DynamicList<scalar> displacementHistory_;

        //- Force sample history
        DynamicList<scalar> forceHistory_;

        //- History file ptr
        autoPtr<OFstream> historyFilePtr_;


    // Private Member Functions

        //- Sample the monitored solid point displacement
        bool sampleDisplacement(scalar& displacement) const;

        //- Sample the monitored interface force
        bool sampleForce(scalar& force) const;

        //- Open the history file
        void openHistoryFile();

        //- Write one history row
        void writeHistory
        (
            const scalar displacement,
            const scalar force,
            const scalar displacementAbsChange,
            const scalar displacementRelChange,
            const scalar forceAbsChange,
            const scalar forceRelChange,
            const bool steady
        );

        //- Execute the steady-state check
        bool checkSteadyState();

        //- Disallow default bitwise copy construct
        fsiSteadyStateControl(const fsiSteadyStateControl&);

        //- Disallow default bitwise assignment
        void operator=(const fsiSteadyStateControl&);


public:

    //- Runtime type information
    TypeName("fsiSteadyStateControl");


    // Constructors

        //- Construct from components
        fsiSteadyStateControl
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
