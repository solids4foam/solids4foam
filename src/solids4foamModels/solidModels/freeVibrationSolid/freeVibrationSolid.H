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
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "solidPolyMesh.H"
#include "BlockSolverPerformance.H"
#include "uniformDimensionedFields.H"

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

        //- Displacement field
        volVectorField D_;

        //- Displacement fields for each mode requested
        PtrList<volVectorField> Dmodes_;

        //- Point mesh
        pointMesh pMesh_;

        //- Point displacement field
        pointVectorField pointD_;

        //- Stress field
        volSymmTensorField sigma_;

        //- Cell-to-point interpolation
        newLeastSquaresVolPointInterpolation volToPoint_;

        //- Cell displacement gradient
        volTensorField gradD_;

        //- Density
        const volScalarField rho_;

        //- Hooke's law second Lame parameter, aka shear modulus
        surfaceScalarField muf_;

        //- Hooke's law first Lame parameter
        surfaceScalarField lambdaf_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField g_;


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

            //- Return const access to the cell displacement field
            const volVectorField& D() const
            {
                return D_;
            }

            //- Return non-const access to thecell displacement field
            volVectorField& D()
            {
                return D_;
            }

            //- Return point mesh
            const pointMesh& pMesh() const
            {
                return pMesh_;
            }

            //- Stress tensor
            const volSymmTensorField& sigma() const
            {
                return sigma_;
            }

            //- Return point displacement field
            virtual const pointVectorField& pointD() const
            {
                return pointD_;
            }

            //- Return multi-material interface pointer
            const newLeastSquaresVolPointInterpolation& volToPoint() const
            {
                return volToPoint_;
            }

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                const label zoneID
            ) const;

            //- Face zone velocity surface (tangential) gradient
            virtual tmp<tensorField> faceZoneSurfaceGradientOfVelocity
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone current normal
            virtual tmp<vectorField> faceZoneNormal
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone primitive patch in current configuration
            virtual tmp<vectorField> currentFaceZonePoints
            (
                const label zoneID
            ) const;

        // Edit

            //- Set traction on the specified patch
            //  This function enables the FSI procedure to pass a traction from
            //  fluid to the solid
            virtual void setTraction
            (
                const label patchID,
                const vectorField& traction
            );

            //- Set pressure at specified patch
            virtual void setPressure
            (
                const label patchID,
                const scalarField& pressure
            );

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& zoneTraction
            );

            //- Set pressure at specified patch
            virtual void setPressure
            (
                const label patchID,
                const label zoneID,
                const scalarField& zonePressure
            );

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            //  Given the user-specified boundary patch traction, this function
            //  should return the surface normal gradient
            //  This function is used by traction-type boundary conditions to
            //  set the boundary gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- End of the program
            virtual void end();

            //- Return const acess to the solution (displacement field)
            virtual const volVectorField& solution() const
            {
                return D_;
            }

            //- Return non-const access to the solution (displacement field)
            virtual volVectorField& solution()
            {
                return D_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
