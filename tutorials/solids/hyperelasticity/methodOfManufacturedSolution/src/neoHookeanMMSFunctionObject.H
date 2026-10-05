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
    Foam::neoHookeanMMSFunctionObject

Description
    Compares the solid model fields with the neo-Hookean manufactured solution
    at the current time. At every time step the mean L1, mean L2 and L-infinity
    norms of the cell displacement error, the point displacement error and the
    cell Cauchy stress error are written to the log. At write times the
    analytical fields and the difference fields are written:

        analyticalD, DDifference
        analyticalPointD, pointDDifference
        analyticalSigma, sigmaDifference

    The solution parameters are read from constant/neoHookeanMMSProperties.

Usage
    \verbatim
    functions
    {
        mms
        {
            type    neoHookeanMMS;
        }
    }
    \endverbatim

SourceFiles
    neoHookeanMMSFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMMSFunctionObject_H
#define neoHookeanMMSFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "neoHookeanMMS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class neoHookeanMMSFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMMSFunctionObject
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Object that calculates the manufactured solution
        autoPtr<neoHookeanMMS> mmsPtr_;


    // Private Member Functions

        //- Write the error norms of a difference field to the log
        template<class Type>
        void writeNorms(const word& title, const Field<Type>& diff) const;

        //- Calculate and write the errors
        bool writeData();

        //- Disallow default bitwise copy construct
        neoHookeanMMSFunctionObject(const neoHookeanMMSFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanMMSFunctionObject&);


public:

    //- Runtime type information
    TypeName("neoHookeanMMS");


    // Constructors

        //- Construct from components
        neoHookeanMMSFunctionObject
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
#else
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
