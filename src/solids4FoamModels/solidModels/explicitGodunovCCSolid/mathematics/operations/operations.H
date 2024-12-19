/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::operations

Description
    Define commonly used mathematical operations.

SourceFiles
    operations.C

\*---------------------------------------------------------------------------*/

#ifndef operations_H
#define operations_H

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class operations
{
    // Private data

        //- Mesh
         const fvMesh& mesh_;

        //- Eigen-values
        vector eigVal_;

        //- Eigen-vectors
        tensor eigVec_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        operations(const operations&);

        //- Disallow default bitwise assignment
        void operator=(const operations&);

public:

    //- Runtime type information
    TypeName("operations");

    // Constructors

        //- Construct from mesh
        operations
        (
            const fvMesh&
        );


    //- Destructor
    virtual ~operations();


    // Member functions

        //Access

            const vector& eigenValue() const
            {
                return eigVal_;
            }

            const tensor& eigenVector() const
            {
                return eigVec_;
            }

        // Edit

            //- Minimum edge length
            dimensionedScalar minimumEdgeLength();

            //- Inverse transpose of a volTensorField
            volTensorField invT
            (
                const GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            //- Tensor cross product of two volTensorFields
            volTensorField tensorProduct
            (
                const GeometricField<tensor, fvPatchField, volMesh>&,
                const GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            //- Tensor cross product of two tensors
            tensor tensorProduct(const tensor&, const tensor&) const;

            //- Decompose a volField tensor
            void decomposeTensor
            (
                const GeometricField<tensor, fvPatchField, volMesh>&,
                GeometricField<vector, fvPatchField, volMesh>&,
                GeometricField<vector, fvPatchField, volMesh>&,
                GeometricField<vector, fvPatchField, volMesh>&
            ) const;

            //- Decompose a surfaceField tensor
            void decomposeTensor
            (
                const GeometricField<tensor, fvsPatchField, surfaceMesh>&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&
            ) const;

            //- Decompose tensor and return x row
            volVectorField decomposeTensorX
            (
                const GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            //- Decompose tensor and return y row
            volVectorField decomposeTensorY
            (
                const GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            //- Decompose tensor and return z row
            volVectorField decomposeTensorZ
            (
                const GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            // Compute eigen-values and eigen-vectors
            void eigenStructure(const tensor&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //