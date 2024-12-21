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

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::mechanics

Description
    Evaluate quantities widely used in continuum mechanics.

SourceFiles
    mechanics.C

\*---------------------------------------------------------------------------*/

#ifndef mechanics_H
#define mechanics_H

#include "fvc.H"
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mechanics
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Operations class
        operations op;

        //- Material normal
        surfaceVectorField N_;

        //- Spatial normal
        surfaceVectorField n_;

        //- Stabilisation matrix for linear momentum
        surfaceTensorField S_lm_;

        //- Stabilisation matrix for traction
        surfaceTensorField S_t_;


        //- CFl number
        scalar cfl_;

        //- Time step
        scalar tStep_;

        //- Stretch
        volScalarField stretch_;

private:

    //- Disallow default bitwise copy construct
        mechanics(const mechanics&);

    //- Disallow default bitwise assignment
        void operator=(const mechanics&);

public:

    //- Runtime type information
    TypeName("mechanics");

    // Constructors

        //- Construct from mesh and dictionary
        mechanics
        (
            const volTensorField&,
            const dictionary&
        );


    //- Destructor
    virtual ~mechanics();


    // Member functions

        // Access

            const surfaceVectorField& spatialNormall() const
            {
                return n_;
            }

            const surfaceTensorField& Smatrix_lm() const
            {
                return S_lm_;
            }

            const surfaceTensorField& Smatrix_t() const
            {
                return S_t_;
            }

            const volScalarField stretch() const
            {
                return stretch_;
            }

        // Edit

            //- Computes spatial surface normals
            surfaceVectorField spatialNormal(const volTensorField&);

            //- Compute stabilisation
            void correct
            (
                const GeometricField<scalar, fvPatchField, volMesh>&,
                const GeometricField<scalar, fvPatchField, volMesh>&,
                const GeometricField<tensor, fvPatchField, volMesh>&
            );

            //-
            void time
            (
                Time&,
                dimensionedScalar&,
                dimensionedScalar
            );

        // Output

            //- Print centroid of geometry
            void printCentroid() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //