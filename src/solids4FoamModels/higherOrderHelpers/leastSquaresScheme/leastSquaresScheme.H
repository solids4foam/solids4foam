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
    Foam::leastSquaresScheme

Description
    Abstract interface for run-time selectable high-order least-squares
    reconstruction schemes.

SourceFiles
    leastSquaresScheme.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresScheme_H
#define leastSquaresScheme_H

#include "runTimeSelectionTables.H"
#include "fvMeshQuadrature.H"
#include "leastSquaresStencil.H"
#include "processorFvPatch.H"
#include "surfaceFields.H"
#include "symmTensor3rdOrder.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class leastSquaresScheme Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresScheme
{
protected:

    // Protected data

        //- Const reference to the mesh
        const fvMesh& mesh_;


    // Protected Member Functions

        //- Calculate value coefficients for one cell at a specified point
        //  Stencil-cell coefficients are followed by the owner coefficient
        virtual void cellValueCoeffsAtPoint
        (
            const label cellI,
            const point& x,
            UList<scalar>& coeffs
        ) const = 0;

        //- Clear cached face-centre value coefficients
        void clearFaceCentreValueCoeffs() const;


private:

    // Private data

        //- Owner-side value coefficients at face centres
        mutable autoPtr<CompactListList<scalar>>
            ownerFaceCentreValueCoeffsPtr_;

        //- Neighbour-side value coefficients at internal-face centres
        mutable autoPtr<CompactListList<scalar>>
            neighbourFaceCentreValueCoeffsPtr_;

    // Private Member Functions

        //- Calculate and cache face-centre value coefficients
        void makeFaceCentreValueCoeffs() const;

        //- For global index return a local or remote field value
        template<class Type, class GlobalCells>
        inline Type fieldValue
        (
            const label globalCellID,
            const GlobalCells& globalCells,
            const Map<FixedList<label, 2>>& remoteLocation,
            const UList<Type>& localField,
            const List<Field<Type>>& remoteField
        ) const;

        //- Apply one owner-last value-coefficient row
        template<class Type, class GlobalCells>
        Type evaluateCellValueCoeffs
        (
            const label cellI,
            const UList<scalar>& coeffs,
            const GlobalCells& globalCells,
            const Map<FixedList<label, 2>>& remoteLocation,
            const UList<Type>& localField,
            const List<Field<Type>>& remoteField
        ) const;

        //- Exchange reconstructed values across a processor patch
        template<class Type>
        void exchangeProcessorPatchField
        (
            const processorFvPatch& procPatch,
            const Field<Type>& sendField,
            Field<Type>& receiveField
        ) const;

        //- Evaluate owner- and neighbour-side values at face centres
        template<class Type>
        void evaluateFaceCentreValues
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            GeometricField<Type, fvsPatchField, surfaceMesh>& ownerValues,
            GeometricField<Type, fvsPatchField, surfaceMesh>& neighbourValues
        ) const;

        //- Disallow default bitwise copy construct
        leastSquaresScheme(const leastSquaresScheme&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresScheme&);


public:

    //- Runtime type information
    TypeName("leastSquaresScheme");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            leastSquaresScheme,
            dictionary,
            (
                const fvMesh& mesh,
                const boolList& includePatchInStencils,
                const dictionary& dict
            ),
            (mesh, includePatchInStencils, dict)
        );


    // Selectors

        //- Create a least-squares reconstruction scheme
        static autoPtr<leastSquaresScheme> New
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );


    // Constructors

        //- Construct from mesh
        explicit leastSquaresScheme(const fvMesh& mesh);


    //- Destructor
    virtual ~leastSquaresScheme();


    // Member Functions

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return polynomial order of the reconstruction
        virtual label polynomialOrder() const = 0;

        //- Return quadrature data
        virtual const fvMeshQuadrature& quadrature() const = 0;

        //- Return stencil data
        virtual const leastSquaresStencil& stencil() const = 0;

        //- Return addressing associated with face-gradient coefficients
        virtual const CompactListList<label>& faceGradStencil() const = 0;

        //- Return face-gradient reconstruction coefficients
        virtual const List<CompactListList<vector>>&
        faceGradCoeffs() const = 0;

        //- Return owner-side face-centre value coefficients
        const CompactListList<scalar>& ownerFaceCentreValueCoeffs() const;

        //- Return neighbour-side face-centre value coefficients
        const CompactListList<scalar>& neighbourFaceCentreValueCoeffs() const;

        //- Clear owned geometry-dependent data and coefficients
        virtual void clear() const = 0;


    // Evaluation

        //- Evaluate a scalar field reconstructed from a cell at a point
        virtual scalar valueAtPoint
        (
            const volScalarField& vf,
            const label cellID,
            const point& x
        ) const = 0;

        //- Evaluate a vector field reconstructed from a cell at a point
        virtual vector valueAtPoint
        (
            const volVectorField& vf,
            const label cellID,
            const point& x
        ) const = 0;

        //- Evaluate owner- and neighbour-side scalar values at face centres
        void faceCentreValues
        (
            const volScalarField& vf,
            surfaceScalarField& ownerValues,
            surfaceScalarField& neighbourValues
        ) const;

        //- Evaluate owner- and neighbour-side vector values at face centres
        void faceCentreValues
        (
            const volVectorField& vf,
            surfaceVectorField& ownerValues,
            surfaceVectorField& neighbourValues
        ) const;

        //- Gradient of a scalar field at face quadrature points
        virtual void fGrad
        (
            const volScalarField& vf,
            CompactListList<vector>& result
        ) const = 0;

        //- Gradient of a vector field at face quadrature points
        virtual void fGrad
        (
            const volVectorField& vf,
            CompactListList<tensor>& result
        ) const = 0;

        //- Scalar field values at patch face quadrature points
        virtual autoPtr<CompactListList<scalar>> patchFaceQuadValues
        (
            const volScalarField& vf,
            const label patchI
        ) const = 0;

        //- Vector field values at patch face quadrature points
        virtual autoPtr<CompactListList<vector>> patchFaceQuadValues
        (
            const volVectorField& vf,
            const label patchI
        ) const = 0;

        //- Cell-centre gradient of a scalar field
        virtual tmp<volVectorField> grad(const volScalarField& vf) const = 0;

        //- Cell-centre gradient of a vector field
        virtual tmp<volTensorField> grad(const volVectorField& vf) const = 0;

        //- Cell-centre second gradient of a scalar field
        virtual tmp<volSymmTensorField> secondGrad
        (
            const volScalarField& vf
        ) const = 0;

        //- Cell-centre third gradient of a scalar field
        virtual autoPtr<List<symmTensor3rdOrder>> thirdGrad
        (
            const volScalarField& vf
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
