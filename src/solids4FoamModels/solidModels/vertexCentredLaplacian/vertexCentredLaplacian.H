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
    vertexCentredLaplacian

Description
    Vertex-centred solver for the Laplace equation.

    Discretisation is performed using a vertex-centred finite volume approach
    where the governing equations are integrated over a dual mesh. The dual mesh
    is automatically constructed from the primary mesh using the meshDualiser.

Author
    Philip Cardiff, UCD.  All rights reserved.
    Federico Mazzanti, UCD. All rights reserved.

SourceFiles
    vertexCentredLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef vertexCentredLaplacian_H
#define vertexCentredLaplacian_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "sparseMatrix.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "globalPointIndices.H"
#ifdef OPENFOAM_COM
    #include "pointVolInterpolation.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                   Class vertexCentredLaplacian Declaration
\*---------------------------------------------------------------------------*/

class vertexCentredLaplacian
:
    public solidModel
{
    // Private Data

        //- Point solution field, e.g. temperature
        pointScalarField pointT_;

        //- Diffusivity
        const scalar diffusivity_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

        //- Field of flags indicating the degrees of freedom where the value is
        //  set
        //  Currently the size of this list is equal to the number of points
        boolList fixedDofs_;

        //- Values for fixed degrees of freedom
        //  The size of this list is equal to the number of points
        scalarField fixedDofValues_;

        //- Scale factor for the fixed degree of freedom equations to improve the
        //  matrix conditioning
        const scalar fixedDofScale_;

        //- Point volume field
        //  Dual mesh cell volumes associated with the primary mesh points
        pointScalarField pointVol_;

        //- Solution gradient field at the dual mesh faces
        surfaceVectorField dualGradTf_;

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;

#ifdef OPENFOAM_COM
        //- Interpolator from points to cells
        pointVolInterpolation pointVolInterp_;
#endif


    // Private Member Functions

        //- Update the source vector for the linear system
        void updateSource
        (
            scalarField& source,
            const pointScalarField& pointT,
            const surfaceVectorField& dualGradTf,
            const scalar diffusivity,
            const labelList& dualCellToPoint
        );

        //- Set fixed degree of freedom lists
        void setFixedDofs
        (
            const pointScalarField& pointT,
            boolList& fixedDofs,
            scalarField& fixedDofValues
        ) const;

        //- Enforce specified gradient boundary conditions on the flux, i.e.
        //  Neumann/natural boundary condition
        void enforceGradientBoundaries
        (
            const pointScalarField& pointT,
            surfaceScalarField& dualFlux,
            const fvMesh& mesh, // corresponding to pointT mesh
            const labelListList& pointToDualFaces
        ) const;

        //- Disallow default bitwise copy construct
        vertexCentredLaplacian(const vertexCentredLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const vertexCentredLaplacian&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            // Not used but required by solidModel base class
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("vertexCentredLaplacian");

    // Constructors

        //- Construct from components
        vertexCentredLaplacian
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~vertexCentredLaplacian();


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                notImplemented("For vertex-centred approach");

                // Keep compiler happy
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
            ) const
            {
                notImplemented("For vertex-centred approach");

                // Keep compiler happy
                return tmp<vectorField>(new vectorField());
            }

        // Write

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
