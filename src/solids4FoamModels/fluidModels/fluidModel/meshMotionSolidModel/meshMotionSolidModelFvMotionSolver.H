/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::meshMotionSolidModelFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh using a solids4foam solidModel.

SourceFiles
    meshMotionSolidModelFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef meshMotionSolidModelFvMotionSolver_H
#define meshMotionSolidModelFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class meshMotionSolidModelFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class meshMotionSolidModelFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Solid model
        mutable autoPtr<solidModel> modelPtr_;


    // Private Member Functions

        //- No copy construct
        meshMotionSolidModelFvMotionSolver
        (
            const meshMotionSolidModelFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const meshMotionSolidModelFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("solids4foamSolidModel");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        meshMotionSolidModelFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        meshMotionSolidModelFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~meshMotionSolidModelFvMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
