/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::totalVelocityLaplacianFvMotionSolver

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion velocity.

SourceFiles
    totalVelocityLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef totalVelocityLaplacianFvMotionSolver_H
#define totalVelocityLaplacianFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
              Class totalVelocityLaplacianFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class totalVelocityLaplacianFvMotionSolver
:
    public fvMotionSolver
{
    // Private data

        //- Point motion field
        mutable pointVectorField pointMotionU_;

        //- Cell-centre motion field
        mutable volVectorField cellMotionU_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        totalVelocityLaplacianFvMotionSolver
        (
            const totalVelocityLaplacianFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const totalVelocityLaplacianFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("totalVelocityLaplacian");


    // Constructors

        //- Construct from polyMesh and data stream
        totalVelocityLaplacianFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    // Destructor

        virtual ~totalVelocityLaplacianFvMotionSolver();


    // Member Functions

        //- Return reference to the point motion velocity field
        pointVectorField& pointMotionU()
        {
            return pointMotionU_;
        }

        //- Return const reference to the point motion velocity field
        const pointVectorField& pointMotionU() const
        {
            return pointMotionU_;
        }

        //- Return reference to the cell motion velocity field
        volVectorField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return const reference to the cell motion velocity field
        const volVectorField& cellMotionU() const
        {
            return cellMotionU_;
        }

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
