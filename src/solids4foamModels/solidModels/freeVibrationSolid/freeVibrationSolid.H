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
    freeVibrationSolid

Description
    Solid model for performing free vibrationn modal analysis.

    This is work in progress.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    freeVibrationSolid.C

\*---------------------------------------------------------------------------*/

#ifndef freeVibrationSolid_H
#define freeVibrationSolid_H

#include "solidModel.H"
#include "solidPolyMesh.H"
#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class freeVibrationSolid Declaration
\*---------------------------------------------------------------------------*/

class freeVibrationSolid
:
    public solidModel
{
    // Private data

        //- The number of models to be written out
        const int nModes_;

        //- The solidPolyMesh stores the additional addressing required for the
        //  extended implicit stencil
        solidPolyMesh extendedMesh_;

        //- Solution vector
        //  This is a field containing the cell displacements followed by the
        //  boundary face displacements
        vectorIOField solutionVec_;

        //- Hooke's law second Lame parameter, aka shear modulus
        surfaceScalarField muf_;

        //- Hooke's law first Lame parameter
        surfaceScalarField lambdaf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        freeVibrationSolid
        (
            const freeVibrationSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const freeVibrationSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("freeVibration");

    // Constructors

        //- Construct from components
        freeVibrationSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~freeVibrationSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                return D();
            }

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
