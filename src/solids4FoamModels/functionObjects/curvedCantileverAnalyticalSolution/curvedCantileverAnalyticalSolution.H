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
    curvedCantileverAnalyticalSolution

Description
    FunctionObject that generates the analytical solution for a curved 
    cantilever beam with end loading
    
    The solution is taken from 
        Sadd MH. Elasticity: Theory, Applications, and Numerics. Elsevier 2009.

    The analytical solution requires the following inputs:
        - rInner : inner beam radius
        - rOuter : outer beam radius
        - force : applied force (in N/m) at beam free end
        - E : Young's modulus
        - nu : Poisson's ratio

    Young's modulus and Poisson's factor are currently excess, but they are left 
    because they will be needed if one wants to add the analytical displacement 
    
    It is assumed that the inner and outer beam surfaces are traction-free.

Author
    Ivan Batistic & Philip Cardiff, UCD. All rights reserved.

SourceFiles
    curvedCantileverAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef curvedCantileverAnalyticalSolution_H
#define curvedCantileverAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class curvedCantileverAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class curvedCantileverAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Inner beam radius
        const scalar rInner_;

        //- Outer beam radius
        const scalar rOuter_;
        
        //- Beam end traction
        const scalar force_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;


    // Private Member Functions
    
        //- Calculate stress tensor for given point
        symmTensor curvedCantileverStress(const vector& C);

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        curvedCantileverAnalyticalSolution
        (
            const curvedCantileverAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const curvedCantileverAnalyticalSolution&);


public:

    //- Runtime type information
    TypeName("curvedCantileverAnalyticalSolution");


    // Constructors

        //- Construct from components
        curvedCantileverAnalyticalSolution
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

#ifdef OPENFOAMESIORFOUNDATION
        //- Write
        virtual bool write();
#endif

#ifndef OPENFOAMESIORFOUNDATION
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
