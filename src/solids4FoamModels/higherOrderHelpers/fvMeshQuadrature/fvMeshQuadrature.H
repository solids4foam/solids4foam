/*---------------------------------------------------------------------------* \
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
    fvMeshQuadrature

Description
    Provide quadrature points and weights for faces and cells. The number of
    quadrature points corresponds to the number of points to exactly integrate
    chosen order of polynomial for cell and face. Points are in global
    coordinates and weights are phsical weights, i.e. they are multiplied with
    triangle area. Becouse of this we don't need to normalise them and multiply
    later with area.

    Face quadrature:
        3D: Cell faces are triangulated using fan triangulation;
            Quadrature points are calculated for each triangle sub-element.
            In case of poor triangulation and active fallback switch algorithm
            will use central point triangulation.

        2D: Mesh is treated as planar; quadrature points are on edges.
            No triangulation needed.

    Cell quadrature:
        3D: Polyhedral cells decomposed into tetrahedra using existing face
            triangulation and cell centroid point.
            Tetrahedral cells are not decomposed.

        2D: Mesh is treated as planar; cell quadrature points are on the plane.
            Planar cells (faces) are triangulated using fan triangulation..

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

SourceFiles
    fvMeshQuadrature.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshQuadrature_H
#define fvMeshQuadrature_H

#include "fvMesh.H"
#include "CompactListList.H"
#include "compatibilityFunctions.H"
#include "symmTensor3rdOrder.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshQuadrature Declaration
\*---------------------------------------------------------------------------*/

class fvMeshQuadrature
{
private:
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Integration orders for cell
        const label cellOrder_;

        //- Integration order for face
        const label faceOrder_;

        //- Minimum area ratio of triangle subelement area and face area
        //  below which decomposition method is changed from fan triangulation to
        //  triangulation using geometric centre
        static const scalar minTriAreaRatio_;

        //-  Minimum ratio of triangle and circum-circle area,
        //   scaled so that an equilateral triangle has a quality of 1
        static const scalar minTriQuality_;

        //- Use triangulation using geometric centre instead of fan
        //  triangulation when minTriAreaRatio is lower than threshold value
        const Switch allowDegenerateTriFallback_;

        //- Face quadrature points
        mutable autoPtr<CompactListList<point>> faceQuadPointsPtr_;

        //- Face quadrature points weights
        mutable autoPtr<CompactListList<scalar>> faceQuadWeightsPtr_;

        //- Cell quadrature points
        mutable autoPtr<CompactListList<point>> cellQuadPointsPtr_;

        //- Cell quadrature points weights
        mutable autoPtr<CompactListList<scalar>> cellQuadWeightsPtr_;

        //- Volume-normalised first-order central cell moments
        mutable autoPtr<List<vector>> firstOrderCellMomentsPtr_;

        //- Volume-normalised second-order central cell moments
        mutable autoPtr<List<symmTensor>> secondOrderCellMomentsPtr_;

        //- Volume-normalised third-order central cell moments
        mutable autoPtr<List<symmTensor3rdOrder>> thirdOrderCellMomentsPtr_;

        //- Disallow default bitwise copy construct
        fvMeshQuadrature(const fvMeshQuadrature&);

        //- Disallow default bitwise assignment
        void operator=(const fvMeshQuadrature&);

        //- Check face triangulation quality, return true if triangulation
        //  have area ratio or quality below treshold values
        bool poorTriangulation
        (
#ifdef OPENFOAM_ORG
            const triFaceList& triFaces,
#else
            const faceList& triFaces,
#endif
            const pointField& pts
        ) const;

        //- Calculate quadrature points and corresponding weights in 3D
        void calcQuadPointsAndWeights3D() const;

        //- Calculate quadrature points and corresponding weights in 2D
        void calcQuadPointsAndWeights2D() const;

        //- Calls the corresponding function for 3D or for 2D meshes
        void calcQuadPointsAndWeights() const;

        //- Calculate volume-normalised first-order cell moments
        void calcFirstOrderCellMoments() const;

        //- Calculate volume-normalised cell moments
        void calcCellMoments() const;

public:
    //- Runtime type information
    TypeName("fvMeshQuadrature");

    // Constructors

        //- Construct from components
        fvMeshQuadrature
        (
            const fvMesh& mesh,
            const label cellOrder,
            const label faceOrder,
            const Switch allowDegenerateTriFallback = true
        );

    // Destructor

        virtual ~fvMeshQuadrature();


    // Member Functions

            //- Return reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return order of accuracy for cell quadrature
            label cellOrder() const
            {
                return cellOrder_;
            }

            //- Return order of accuracy for face quadrature
            label faceOrder() const
            {
                return faceOrder_;
            }

            //- Return the face quadrature weights
            const CompactListList<scalar>& faceQuadWeights() const;

            //- Return the face quadrature points
            const CompactListList<point>& faceQuadPoints() const;

            //- Return the cell quadrature weights
            const CompactListList<scalar>& cellQuadWeights() const;

            //- Return the cell quadrature points
            const CompactListList<point>& cellQuadPoints() const;

            //- Return volume-normalised first-order cell moments
            const List<vector>& firstOrderCellMoments() const;

            //- Return volume-normalised second-order cell moments
            const List<symmTensor>& secondOrderCellMoments() const;

            //- Return volume-normalised third-order cell moments
            const List<symmTensor3rdOrder>& thirdOrderCellMoments() const;

            //- Clear data after topo change
            void clear();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
