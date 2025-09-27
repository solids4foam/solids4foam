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
    Foam::decayingTaylorGreenVortexMeshMotion

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver which prescribes the mesh motion for the
    decayingTaylorGreenVortex tutorial case. The purpose of this solver is for
    testing the performance of a moving mesh fluid solver. The mesh motion is
    given by

        x(x,y,t) = y(x,y,t) = A*sin(π*t/0.8)*sin(π*x)*sin(π*y)

    where the scaling factor A defaults to 0.1 but can be supplied by the user.

SourceFiles
    decayingTaylorGreenVortexMeshMotion.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef decayingTaylorGreenVortexMeshMotion_H
#define decayingTaylorGreenVortexMeshMotion_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class decayingTaylorGreenVortexMeshMotion Declaration
\*---------------------------------------------------------------------------*/

class decayingTaylorGreenVortexMeshMotion
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Motion scale factor
        const scalar A_;

        //- Flag to indicate orthogonal mesh motion or more general motion
        const Switch orthogonalMeshMotion_;


    // Private Member Functions

        //- No copy construct
        decayingTaylorGreenVortexMeshMotion
        (
            const decayingTaylorGreenVortexMeshMotion&
        ) = delete;

        //- No copy assignment
        void operator=(const decayingTaylorGreenVortexMeshMotion&) = delete;


public:

    //- Runtime type information
    TypeName("decayingTaylorGreenVortexMeshMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        decayingTaylorGreenVortexMeshMotion
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        decayingTaylorGreenVortexMeshMotion
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~decayingTaylorGreenVortexMeshMotion();


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
