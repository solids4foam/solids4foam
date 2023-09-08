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
    contactPatchTestAnalyticalSolution

Description
    FunctionObject that generates the analytical solution for contact patch test

    The solution is taken from
        Crisfield MA. Re-visiting the contact patch test.
        Int J Numer Methods Eng. 2000

    The analytical solution requires the following inputs:
        - displacement : upper block top surface vertical displacement
        - E : Young's modulus
        - nu : Poisson's ratio

Author
    Ivan Batistic & Philip Cardiff, UCD. All rights reserved.

SourceFiles
    contactPatchTestAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef contactPatchTestAnalyticalSolution_H
#define contactPatchTestAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class contactPatchTestAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class contactPatchTestAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Beam end traction
        const scalar displacement_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;


    // Private Member Functions

        //- Calculate stress tensor for given point
        symmTensor contactPatchTestStress();

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        contactPatchTestAnalyticalSolution
        (
            const contactPatchTestAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const contactPatchTestAnalyticalSolution&);


public:

    //- Runtime type information
    TypeName("contactPatchTestAnalyticalSolution");


    // Constructors

        //- Construct from components
        contactPatchTestAnalyticalSolution
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
