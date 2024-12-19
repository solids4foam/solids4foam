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
    Foam::gradientSchemes

Description
    Define gradient computation schemes.

SourceFiles
    gradientSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gradientSchemes_H
#define gradientSchemes_H

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class gradientSchemes
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Cell owners
        const labelUList& own_;

        //- Cell neighbours
        const labelUList& nei_;

        // Cell centre coordinates
        const volVectorField& X_;

        // Face centre coordinates
        const surfaceVectorField& XF_;

        //- Material nodal coordinates
        const pointField& XN_;

        //- Inverse of distance matrix
        volTensorField Ainv_;

        //- Inverse of local distance matrix
        volTensorField AinvLocal_;

private:

    // Private member functions

        //- Disallow default bitwise copy construct
        gradientSchemes(const gradientSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const gradientSchemes&);


public:

    //- Runtime type information
    TypeName("gradientSchemes");

    // Constructors

        //- Construct from mesh
        gradientSchemes(const fvMesh&);


    // Destructor
    virtual ~gradientSchemes();


    // Member functions

        // Edit

            //- Inverse distance matrix for gradient
            void distanceMatrix
            (
                GeometricField<tensor, fvPatchField, volMesh>&
            );

            //- Inverse distance matrix for local gradient
            void distanceMatrixLocal
            (
                GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            //- Least square gradient for a volScalarField
            volVectorField gradient
            (
                const GeometricField<scalar, fvPatchField, volMesh>& U
            ) const;

            //- Least square gradient for a volVectorField
            volTensorField gradient
            (
                const GeometricField<vector, fvPatchField, volMesh>& U
            ) const;

            //- Least square gradient for a volTensorField
            void gradient
            (
                const GeometricField<tensor, fvPatchField, volMesh>&,
                GeometricField<tensor, fvPatchField, volMesh>&,
                GeometricField<tensor, fvPatchField, volMesh>&,
                GeometricField<tensor, fvPatchField, volMesh>&
            ) const;

            //- Local least square gradient for a volVectorField
            volTensorField localGradient
            (
                const GeometricField<vector, fvPatchField, volMesh>&,
                const GeometricField<vector, fvsPatchField, surfaceMesh>&
            ) const;

            //- Reconstruction of a volScalarField to surfaceScalarField
            void reconstruct
            (
                GeometricField<scalar, fvPatchField, volMesh>&,
                const GeometricField<vector, fvPatchField, volMesh>&,
                GeometricField<scalar, fvsPatchField, surfaceMesh>&,
                GeometricField<scalar, fvsPatchField, surfaceMesh>&
            );

            //- Reconstruction of a volVectorField to surfaceVectorField
            void reconstruct
            (
                GeometricField<vector, fvPatchField, volMesh>&,
                const GeometricField<tensor, fvPatchField, volMesh>&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&
            );

            //- Reconstruction for a volTensorField to surfaceTensorField
            void reconstruct
            (
                GeometricField<tensor, fvPatchField, volMesh>&,
                const GeometricField<tensor, fvPatchField, volMesh>&,
                const GeometricField<tensor, fvPatchField, volMesh>&,
                const GeometricField<tensor, fvPatchField, volMesh>&,
                GeometricField<tensor, fvsPatchField, surfaceMesh>&,
                GeometricField<tensor, fvsPatchField, surfaceMesh>&
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //