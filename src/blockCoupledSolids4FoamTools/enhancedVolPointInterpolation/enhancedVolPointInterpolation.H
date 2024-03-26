/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::enhancedVolPointInterpolation

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting and the gradient extrapolated values.

    The value at the points is extrapolated from the surrounding cell centres
    using the cell-centre values and the cell-centre gradients. The inverse
    distance weighted average is then taken at the point.

    This approach was proposed in:
    Z. Tukovi c and H. Jasak. Updated lagrangian finite volume solver for large
    deformation dynamic response of elastic body. Transactions of FAMENA,
    1(31):1–16, 2007.

SourceFiles
    enhancedVolPointInterpolation.C
    enhancedVolPointInterpolate.C

Author
    Extension to gradient approach by Philip Cardiff based on the method above

\*---------------------------------------------------------------------------*/

#ifndef enhancedVolPointInterpolation_H
#define enhancedVolPointInterpolation_H

#include "MeshObject.H"
#include "scalarList.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class enhancedVolPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class enhancedVolPointInterpolation
:
    public MeshObject<fvMesh, UpdateableMeshObject, enhancedVolPointInterpolation>
{
    // Private data

        //- Interpolation scheme weighting factor array.
        scalarListList pointWeights_;


    // Boundary handling

        //- Boundary addressing
        autoPtr<primitivePatch> boundaryPtr_;

        //- Per boundary face whether is on non-coupled, non-empty patch
        bitSet boundaryIsPatchFace_;

        //- Per mesh(!) point whether is on non-coupled, non-empty patch (on
        //  any processor)
        bitSet isPatchPoint_;

        //- Per boundary point the weights per pointFaces.
        scalarListList boundaryPointWeights_;

        //- Whether mesh has any non-processor coupled patches
        bool hasSeparated_;

        //- Optional scaling for coupled point patch fields using the
        //  swapAddSeparated to add additional contributions (e.g. cyclicAMI).
        //  (note: one-to-one matching (processor) uses a different mechanism)
        tmp<scalarField> normalisationPtr_;


    // Private Member Functions

        //- Constructor helper: check if any non-processor coupled patches
        static bool hasSeparated(const pointMesh& pMesh);

        //- Construct addressing over all boundary faces
        void calcBoundaryAddressing();

        //- Make weights for internal and coupled-only boundarypoints
        void makeInternalWeights(scalarField& sumWeights);

        //- Make weights for points on uncoupled patches
        void makeBoundaryWeights(scalarField& sumWeights);

        //- Construct all point weighting factors
        void makeWeights();

        //- Helper: interpolate oneField onto boundary points only. Used to
        //  correct normalisation
        void interpolateOne
        (
            const tmp<scalarField>& tnormalisation,
            pointScalarField& pf
        ) const;

        //- Helper: push master point data to collocated points
        template<class Type>
        void pushUntransformedData(List<Type>&) const;

        //- Get boundary field in same order as boundary faces. Field is
        //  zero on all coupled and empty patches
        template<class Type>
        tmp<Field<Type>> flatBoundaryField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Add separated contributions
        template<class Type>
        void addSeparated
        (
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- No copy construct
        enhancedVolPointInterpolation(const enhancedVolPointInterpolation&) = delete;

        //- No copy assignment
        void operator=(const enhancedVolPointInterpolation&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("enhancedVolPointInterpolation");


    // Constructors

        //- Constructor given fvMesh and pointMesh.
        explicit enhancedVolPointInterpolation(const fvMesh&);


    //- Destructor
    ~enhancedVolPointInterpolation();


    // Member functions

    // Edit

        //- Update mesh topology using the morph engine
        void updateMesh(const mapPolyMesh&);

        //- Correct weighting factors for moving mesh.
        bool movePoints();


    // Interpolation Functions

        //- Interpolate from volField to pointField
        template<class Type>
        void interpolate
        (
            const DimensionedField<Type, volMesh>& vf,
            const DimensionedField
            <
                typename outerProduct<vector, Type>::type, volMesh
            >& gradVf,
            DimensionedField<Type, pointMesh>& pf
        ) const;

        //- Interpolate volField using inverse distance weighting
        //  \return pointField
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate tmp<volField> using inverse distance weighting
        //  \return pointField
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;

        //- Interpolate volField using inverse distance weighting
        //  returning pointField with the same patchField types. Assigns
        //  to any fixedValue boundary conditions to make them consistent
        //  with internal field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const wordList& patchFieldTypes
        ) const;

        //- Interpolate tmp<volField> using inverse distance weighting
        //  returning pointField with the same patchField types. Assigns
        //  to any fixedValue boundary conditions to make them consistent
        //  with internal field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
            const wordList& patchFieldTypes
        ) const;

        //- Interpolate dimensionedField using inverse distance weighting
        //  \return pointField
        template<class Type>
        tmp<DimensionedField<Type, pointMesh>> interpolate
        (
            const DimensionedField<Type, volMesh>&
        ) const;

        //- Interpolate tmp<dimensionedField> using inverse distance
        //  weighting returning pointField
        template<class Type>
        tmp<DimensionedField<Type, pointMesh>> interpolate
        (
            const tmp<DimensionedField<Type, volMesh>>&
        ) const;


    // Low level

        //- Interpolate internal field from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void interpolateInternalField
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate boundary field without applying constraints/boundary
        //  conditions
        template<class Type>
        void interpolateBoundaryField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,

            GeometricField<Type, pointPatchField, pointMesh>& pf
        ) const;

        //- Interpolate boundary with constraints/boundary conditions
        template<class Type>
        void interpolateBoundaryField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            GeometricField<Type, pointPatchField, pointMesh>& pf,
            const bool overrideFixedValue
        ) const;

        //- Interpolate from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate volField using inverse distance weighting
        //  returning pointField with name. Optionally caches
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const word& name,
            const bool cache
        ) const;

        //- Interpolate dimensioned internal field from cells to points
        //  using inverse distance weighting
        template<class Type>
        void interpolateDimensionedInternalField
        (
            const DimensionedField<Type, volMesh>& vf,
            DimensionedField<Type, pointMesh>& pf
        ) const;

        //- Interpolate dimensionedField using inverse distance weighting
        //  returning pointField with name. Optionally caches
        template<class Type>
        tmp<DimensionedField<Type, pointMesh>> interpolate
        (
            const DimensionedField<Type, volMesh>&,
            const word& name,
            const bool cache
        ) const;

    // Interpolation for displacement (applies 2D correction)

        //- Interpolate from volField to pointField
        //  using inverse distance weighting
        void interpolateDisplacement
        (
            const volVectorField&,
            pointVectorField&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "enhancedVolPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
