/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    cantileverAnalyticalSolution

Description
    FunctionObject that generates the analytical fields for the bending slender
    cantilever problem.

    The solution is taken from:
    C.E. Augarde, A.J. Deeks, The use of Timoshenko’s exact solution for a
    cantilever beam in adaptive analysis. Finite Elements in Analysis and
    Design, 44, 2008, 595–601, 10.1016/j.finel.2008.01.010.

    It is assumed that the origin is at the centre of the fixed patch.

    The beam is L long and D deep (thick).

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    cantileverAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef cantileverAnalyticalSolution_H
#define cantileverAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cantileverAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class cantileverAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Load applied in the minus y direction at the other end of the beam
        const scalar P_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Length of the beam
        const scalar L_;

        //- Depth of the beam
        const scalar D_;

        //- Second moment of area
        const scalar I_;

    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        cantileverAnalyticalSolution
        (
            const cantileverAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const cantileverAnalyticalSolution&);


public:

    //- Runtime type information
    TypeName("cantileverAnalyticalSolution");


    // Constructors

        //- Construct from components
        cantileverAnalyticalSolution
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
