/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    coupledUnsLinGeomLinearElasticSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed.

    A block-coupled solution methodology is employed, in contrast to the
    segregated solution methodology employed in linGeomSolid.

    Only Hookean linear elastic material law is allowed; this solidModel
    corresponds closely to the approach used in the reference:
    P. Cardiff, Z Tukovic, H. Jasak, and A. Ivanković. A Block-Coupled Finite
    Volume Methodology for Linear Elasticity and Unstructured Meshes, Computers
    and Structures, 2016, 10.1016/j.compstruc.2016.07.004.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    coupledUnsLinGeomLinearElasticSolid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledUnsLinGeomLinearElasticSolid_H
#define coupledUnsLinGeomLinearElasticSolid_H

#include "solidModel.H"
#ifdef FOAMEXTEND
    #include "solidPolyMesh.H"
    #include "BlockSolverPerformance.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
               Class coupledUnsLinGeomLinearElasticSolid Declaration
\*---------------------------------------------------------------------------*/

class coupledUnsLinGeomLinearElasticSolid
:
    public solidModel
{
    // Private data

#ifdef FOAMEXTEND
        //- The solidPolyMesh stores the additional addressing required for the
        //  extended implicit stencil
        solidPolyMesh extendedMesh_;

        //- Solution vector
        //  This is a field containing the cell displacements followed by the
        //  boundary face displacements
        vectorIOField solutionVec_;
#endif

        //- Hooke's law second Lame parameter, aka shear modulus
        surfaceScalarField muf_;

        //- Hooke's law first Lame parameter
        surfaceScalarField lambdaf_;

    // Private Member Functions

        //- Write banner
        void writeBanner(messageStream& os) const;

        //- Disallow default bitwise copy construct
        coupledUnsLinGeomLinearElasticSolid
        (
            const coupledUnsLinGeomLinearElasticSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const coupledUnsLinGeomLinearElasticSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("coupledUnsLinearGeometryLinearElastic");

    // Constructors

        //- Construct from components
        coupledUnsLinGeomLinearElasticSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );



    // Destructor

        virtual ~coupledUnsLinGeomLinearElasticSolid()
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
