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
    movingLeastSquares

Description
    The class constructs local polynomial approximations of a scalar field
    using a weighted least-squares procedure over face or cell stencils.
    For each reconstruction point (cell centre or face quadrature point),
    a local polynomial expansion of order N is fitted to neighbouring
    stencil cells.

    The resulting coefficients allow reconstruction of field value, gradient,
    second and third derivative. Reconstruction coefficients are precomputed
    and stored for reuse during the simulation, allowing  evaluation of
    derivatives during discretisation phase.

    For each cell centre and face quadrature point the class solves a small
    dense weighted least-squares system of the form A w = b where the matrix A
    contains weighted polynomial basis functions  evaluated at stencil cell
    centres.

    Source:
        Castrillo P, Canelas A, Schillaci E, Rigola J, Oliva A. High-order
        finite volume method for linear elasticity on unstructured meshes.
        Comput Struct, 2022.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

SourceFiles
    movingLeastSquares.C
    movingLeastSquaresTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef movingLeastSquares_H
#define movingLeastSquares_H

#include "leastSquaresScheme.H"
#include "fvMesh.H"
#include "CompactListList.H"
#include "globalIndex.H"
#include <Eigen/Dense>
#include "tensorList.H"
#include "symmTensor3rdOrder.H"
#include "weightFunction.H"
#include "fvMeshQuadrature.H"
#include "compatibilityFunctions.H"
#include "leastSquaresStencil.H"
#include "volFields.H"

#include "fixedDisplacementFvPatchVectorField.H"
#include "processorPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "emptyPolyPatch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class movingLeastSquares Declaration
\*---------------------------------------------------------------------------*/

class movingLeastSquares
:
    public leastSquaresScheme
{
private:
    // Private data

        //- Cell and face stencil pointer
        mutable autoPtr<leastSquaresStencil> stencilPtr_;

        //- Quadrature points and weights pointer
        mutable autoPtr<fvMeshQuadrature> quadraturePtr_;

        //- Weight function pointer
        mutable autoPtr<weightFunction> weightFuncPtr_;

        //- Polynomial order used for the local reconstruction
        const label polynomialOrder_;

        //- Additional cells in the face stencil
        const label faceStencilExtraCells_;

        //- Additional cells in the cell stencil
        const label cellStencilExtraCells_;

        //- Halo depth scaling used when building the stencil
        const scalar haloDepthScale_;

        //- Weight function coefficients dictionary
        const dictionary weightFunctionCoeffs_;

        //- Include patch in the stencils as ghost point
        //  Used for patches with prescribed value
        const boolList includePatchInStencils_;

        //- Calculate the condition number field for cells and faces
        //  Very expensive but useful, use only for debugging purposes
        const Switch calcConditionNumber_;

        //- Cell centre interpolation coefficients
        mutable autoPtr<CompactListList<scalar>>  cellInterpCoeffsPtr_;

        //- Cell centre gradient coefficients
        mutable autoPtr<CompactListList<vector>> cellGradCoeffsPtr_;

        //- Cell centre second derivative tensor coefficients
        mutable autoPtr<CompactListList<symmTensor>> cellSecondGradCoeffsPtr_;

        //- Cell centre third derivative tensor coefficients
        mutable autoPtr
        <
            CompactListList<symmTensor3rdOrder>
        > cellThirdGradCoeffsPtr_;

        //- QR Gradient face Gauss points coefficients
        mutable autoPtr<List<CompactListList<vector>>> faceGradCoeffsPtr_;

        //- Cell condition number field of the least squares system
        mutable autoPtr<volScalarField> cellConditionNumberPtr_;

        //- Face condition number field of the least squares system
        mutable autoPtr<surfaceScalarField> faceConditionNumberPtr_;

    // Helpers

        //- Resulting matrix from QR decomposition (Householder Method)
        struct QRSolution
        {
            Eigen::MatrixXd A;
            scalar cond = 0.0;
        };

        struct derivativeRows
        {
            label ixx, iyy, ixy, izz, ixz, iyz;
            label ixxx, ixxy, ixxz, ixyy, ixyz, ixzz, iyyy, iyyz, iyzz, izzz;

            derivativeRows()
            :
            ixx(-1),
            iyy(-1),
            ixy(-1),
            izz(-1),
            ixz(-1),
            iyz(-1),
            ixxx(-1),
            ixxy(-1),
            ixxz(-1),
            ixyy(-1),
            ixyz(-1),
            ixzz(-1),
            iyyy(-1),
            iyyz(-1),
            iyzz(-1),
            izzz(-1)
            {}
        };

    // Private Member Functions

        //- Return minimum number of points in stencil for polynomialOrder_
        label minNn() const;

        //- Calculate exponents for the Taylor series expansion terms
        void generateExponents
        (
            const label N,
            DynamicList<FixedList<label, 3>>& exponents
        ) const;

        //- Return the row index for exponent (i,j,k) in the generated
        //  Taylor basis
        label rowOf
        (
            const List<FixedList<label, 3>>& exponents,
            label i,
            label j,
            label k
        ) const;

        derivativeRows calcDerivativeRows
        (
            const DynamicList<FixedList<label, 3>>& exponents,
            const bool twoD
        ) const;

         //- Return the cell condition number field
        volScalarField& cellConditionNumber() const;

        //- Return the face condition number field
        surfaceScalarField& faceConditionNumber() const;

        //- Make the cell condition number field
        void makeCellConditionNumber() const;

        //- Make the face condition number field
        void makeFaceConditionNumber() const;

        //- Make the cell and face stencil pointer
        void makeStencils() const;

        //- Make the fvMeshQuadrature pointer
        void makeQuadrature() const;

        //- Make the weight function
        void makeWeightFunction() const;

        //- Calculate interpolation coefficients for the cell centres
        void calcCellCoeffs() const;

        //- Calculate interpolation coefficients for the faces quadrature points
        void calcFaceCoeffs() const;

        //- Calculate value coefficients for one cell at a specified point
        //  Stencil-cell coefficients are followed by the owner coefficient
        virtual void cellValueCoeffsAtPoint
        (
            const label cellI,
            const point& x,
            UList<scalar>& coeffs
        ) const;

        //- Solve using QR decomposition (Householder Method)
        QRSolution QRSolve
        (
            const Eigen::MatrixXd& Q,
            const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& W
        ) const;

        //- Return the cell centre interpolation coefficients
        const CompactListList<scalar>& cellInterpCoeffs() const;

        //- Return the cell centre gradient coefficients
        const CompactListList<vector>& cellGradCoeffs() const;

        //- Return the cell centre second derivative coefficients
        const CompactListList<symmTensor>&  cellSecondGradCoeffs() const;

        //- Return the cell centre  third derivative coefficients
        const CompactListList<symmTensor3rdOrder>& cellThirdGradCoeffs() const;

        //- For global index return field value. Index can refer either to local
        //  or remote field.
        template<class Type, class GlobalCells>
        inline Type fieldValue
        (
            const label globalCellID,
            const GlobalCells& globalCells,
            const Map<FixedList<label, 2>>& remoteLocation,
            const UList<Type>& localField,
            const List<Field<Type>>& remoteField
        ) const;

        //- Evaluate a field using the reconstruction associated with one cell
        template<class Type>
        Type evaluateAtPoint
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const label cellID,
            const point& x
        ) const;

        //- Disallow default bitwise copy construct
        movingLeastSquares(const movingLeastSquares&);

        //- Disallow default bitwise assignment
        void operator=(const movingLeastSquares&);

    public:

    //- Runtime type information
    TypeName("movingLeastSquares");

    // Constructors

        //- Construct from components
        movingLeastSquares
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );

    //- Destructor
    virtual ~movingLeastSquares();


    // Member Functions

        //- Return polynomial order of the local reconstruction
        virtual label polynomialOrder() const
        {
            return polynomialOrder_;
        };

        //- Const access to the weight function
        const weightFunction& weightFunc() const
        {
            if (weightFuncPtr_.empty())
            {
                makeWeightFunction();
            }

            return *weightFuncPtr_;
        }

        //- Const access to the quadrature data
        virtual const fvMeshQuadrature& quadrature() const
        {
            if (quadraturePtr_.empty())
            {
                makeQuadrature();
            }

            return *quadraturePtr_;
        }

        //- Return the cell and face stencil
        virtual const leastSquaresStencil& stencil() const;

        //- Return addressing associated with face-gradient coefficients
        virtual const CompactListList<label>& faceGradStencil() const;

        //- Return the face quadrature points gradient coefficients
        virtual const List<CompactListList<vector>>& faceGradCoeffs() const;

        //- Clear owned MLS data and cached coefficient data
        virtual void clear() const;

    // Evaluation

        //- Evaluate a scalar field reconstructed from a cell at a point
        virtual scalar valueAtPoint
        (
            const volScalarField& vf,
            const label cellID,
            const point& x
        ) const;

        //- Evaluate a vector field reconstructed from a cell at a point
        virtual vector valueAtPoint
        (
            const volVectorField& vf,
            const label cellID,
            const point& x
        ) const;

        //- Gradient of field at face quadrature points
        template<class Type>
        autoPtr<CompactListList<typename outerProduct<vector, Type>::type>>
        fGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Gradient of field at face quadrature points
        template<class Type>
        void fGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            CompactListList<typename outerProduct<vector, Type>::type>& result
        ) const;

        //- Gradient of scalar field at face quadrature points
        virtual void fGrad
        (
            const volScalarField& vf,
            CompactListList<vector>& result
        ) const;

        //- Gradient of vector field at face quadrature points
        virtual void fGrad
        (
            const volVectorField& vf,
            CompactListList<tensor>& result
        ) const;

        //- Get patch face quadrature values
        template<class Type>
        autoPtr<CompactListList<Type>> patchFaceQuadValues
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const label patchI
        ) const;

        //- Get patch face quadrature values for scalar field
        virtual autoPtr<CompactListList<scalar>> patchFaceQuadValues
        (
            const volScalarField& vf,
            const label patchI
        ) const;

        //- Get patch face quadrature values for vector field
        virtual autoPtr<CompactListList<vector>> patchFaceQuadValues
        (
            const GeometricField<vector, fvPatchField, volMesh>& vf,
            const label patchI
        ) const;

        //- Cell centre gradient field
        template<class Type>
        tmp<GeometricField<
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >> grad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Cell centre gradient of scalar field
        virtual tmp<volVectorField> grad(const volScalarField& vf) const;

        //- Cell centre gradient of vector field
        virtual tmp<volTensorField> grad(const volVectorField& vf) const;

        //- Calculate the second grad at the cell centres for the given scalar
        //  field. For vector fields we apply this per component.
        virtual tmp<volSymmTensorField> secondGrad
        (
            const volScalarField& s
        ) const;

        //- Calculate the third grad at the cell centres for the given scalar
        //  field. For vector fields we apply this per component.
        virtual autoPtr<List<symmTensor3rdOrder>>
        thirdGrad(const volScalarField& s) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "movingLeastSquaresTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
