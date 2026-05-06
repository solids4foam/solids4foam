/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    solids4foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    pressurisedCylinderAnalyticalSolution

Description
    FunctionObject that generates analytical displacement and stress fields for
    the internally pressurised cylinder in the cylindrical coordinate system.
    The numerical displacement and stress fields are transformed to cylindrical
    coordinates before the difference fields and norms are calculated.

    Analytical solution taken from:
    J. G. S.P. Timoshenko, Theory of Elasticity. McGraw-Hill, New York, 1970.

    It is assumed that the origin is at (0,0,0). Empty direction is z-axis.

SourceFiles
    pressurisedCylinderAnalyticalSolution.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pressurisedCylinderAnalyticalSolution_H
#define pressurisedCylinderAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "pressurisedCylinderAnalytical.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class pressurisedCylinderAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class pressurisedCylinderAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

       //- Reference to main object registry
        const Time& time_;

        //- Helper object for calculating the analytical displacements and
        //  stress
        pressurisedCylinderAnalytical analyticalSol_;

        //- Calculate cell displacement fields
        const Switch cellDisplacement_;

        //- Calculate point displacement fields
        const Switch pointDisplacement_;

        //- Calculate cell stress fields
        const Switch cellStress_;

        //- Calculate point stress fields
        const Switch pointStress_;

    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        pressurisedCylinderAnalyticalSolution
        (
            const pressurisedCylinderAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const pressurisedCylinderAnalyticalSolution&);

public:

    //- Runtime type information
    TypeName("pressurisedCylinderAnalyticalSolution");


    // Constructors

        //- Construct from components
        pressurisedCylinderAnalyticalSolution
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
