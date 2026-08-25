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
    Foam::kExactLeastSquares

Description
    Run-time selectable skeleton for cell-average k-exact least-squares
    reconstruction. Numerical reconstruction operations are not implemented.

    Sources:
        Sejekan CB, Ollivier-Gooch CF. Improving finite-volume diffusive
        fluxes through better reconstruction. Computers & Fluids, 139:
        216-232, 2016.

        Jalali A, Ollivier-Gooch CF. Higher-Order Finite Volume Solution
        Reconstruction on Highly Anisotropic Meshes. 21st AIAA Computational
        Fluid Dynamics Conference, 2013.

        Pereira da Costa RD. Advanced computational methods towards
        high-performance polymer processing simulations. PhD thesis, 2020.

        Ollivier-Gooch C, Van Altena M. A High-Order-Accurate Unstructured
        Mesh Finite-Volume Scheme for the Advection-Diffusion Equation.
        Journal of Computational Physics, 181:729-752, 2002.

SourceFiles
    kExactLeastSquares.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef kExactLeastSquares_H
#define kExactLeastSquares_H

#include "leastSquaresScheme.H"
#include "processorFvPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "weightFunction.H"
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class kExactLeastSquares Declaration
\*---------------------------------------------------------------------------*/

class kExactLeastSquares
:
    public leastSquaresScheme
{
    // Private data

        //- Cell and face stencil pointer
        mutable autoPtr<leastSquaresStencil> stencilPtr_;

        //- Quadrature points and weights pointer
        mutable autoPtr<fvMeshQuadrature> quadraturePtr_;

        //- Weight function pointer
        mutable autoPtr<weightFunction> weightFuncPtr_;

        //- Polynomial order used for the local reconstruction
        const label polynomialOrder_;

        //- Additional cells in the cell stencil
        const label cellStencilExtraCells_;

        //- Halo depth scaling used when building the stencil
        const scalar haloDepthScale_;

        //- Weight function coefficients dictionary
        const dictionary weightFunctionCoeffs_;

        //- Include patch values in the boundary reconstruction
        const boolList includePatchInStencils_;

        //- Calculate the condition number field for cells
        //  Very expensive but useful, use only for debugging purposes
        const Switch calcConditionNumber_;

        //- Cell centre gradient coefficients
        mutable autoPtr<CompactListList<vector>> cellGradCoeffsPtr_;

        //- Cell centre second derivative tensor coefficients
        mutable autoPtr<CompactListList<symmTensor>> cellSecondGradCoeffsPtr_;

        //- Cell centre third derivative tensor coefficients
        mutable autoPtr
        <
            CompactListList<symmTensor3rdOrder>
        > cellThirdGradCoeffsPtr_;

        //- Addressing associated with face-gradient coefficients
        mutable autoPtr<CompactListList<label>> faceGradStencilPtr_;

        //- Face-gradient reconstruction coefficients
        mutable autoPtr<List<CompactListList<vector>>> faceGradCoeffsPtr_;

        //- Prescribed boundary quadrature-point addressing used by each face
        mutable autoPtr
        <
            List<List<FixedList<label, 2>>>
        > faceBoundaryDataStencilPtr_;

        //- Face-gradient coefficients multiplying prescribed boundary data
        mutable autoPtr
        <
            List<CompactListList<vector>>
        > faceBoundaryDataCoeffsPtr_;

        //- Cell condition number field of the least squares system
        mutable autoPtr<volScalarField> cellConditionNumberPtr_;

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

        //- Return number of non-constant polynomial terms for polynomialOrder_
        label minNn() const;

        //- Calculate exponents for the non-constant Taylor expansion terms
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

        //- Make the cell condition number field
        void makeCellConditionNumber() const;

        //- Make the cell and face stencil pointer
        void makeStencils() const;

        //- Make the fvMeshQuadrature pointer
        void makeQuadrature() const;

        //- Make the weight function
        void makeWeightFunction() const;

        //- Calculate interpolation coefficients for the cell centres
        void calcCellCoeffs() const;

        //- Calculate the cell reconstruction gradient coefficients at a
        //  specified point, aligned with the cell stencil
        void cellGradCoeffsAtPoint
        (
            const label cellI,
            const point& x,
            UList<vector>& coeffs
        ) const;

        //- Calculate value coefficients for one cell at a specified point
        //  Stencil-cell coefficients are followed by the owner coefficient
        virtual void cellValueCoeffsAtPoint
        (
            const label cellI,
            const point& x,
            UList<scalar>& coeffs
        ) const;

        //- Make merged owner- and neighbour-cell stencils for internal faces
        void makeFaceGradStencil() const;

        //- Calculate internal-face gradient coefficients at quadrature points
        //  by merging contribution from owner and neighbour side, and
        //  weighted coefficients for fixed-value boundary faces
        void calcFaceGradCoeffs() const;

        //- Return prescribed boundary quadrature-point addressing
        const List<List<FixedList<label, 2>>>&
        faceBoundaryDataStencil() const;

        //- Return face-gradient coefficients for prescribed boundary data
        const List<CompactListList<vector>>&
        faceBoundaryDataCoeffs() const;

        //- Solve using QR decomposition (Householder Method)
        QRSolution QRSolve
        (
            const Eigen::MatrixXd& Q,
            const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& W
        ) const;

        //- Return the cell centre gradient coefficients
        const CompactListList<vector>& cellGradCoeffs() const;

        //- Return the cell centre second derivative coefficients
        const CompactListList<symmTensor>& cellSecondGradCoeffs() const;

        //- Return the cell centre third derivative coefficients
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

        //- Contract second-derivative and second-moment tensors
        static scalar secondDerivativeMomentContraction
        (
            const symmTensor& secondDerivative,
            const symmTensor& secondMoment
        );

        //- Contract third-derivative and third-moment tensors
        static scalar thirdDerivativeMomentContraction
        (
            const symmTensor3rdOrder& thirdDerivative,
            const symmTensor3rdOrder& thirdMoment
        );

        //- Evaluate a field using the reconstruction associated with one cell
        template<class Type>
        Type evaluateAtPoint
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const label cellID,
            const point& x
        ) const;

        //- Exchange a field across a processor patch
        template<class Type>
        void exchangeProcessorPatchField
        (
            const processorFvPatch& procPatch,
            const Field<Type>& sendField,
            Field<Type>& receiveField
        ) const;

        //- Evaluate the field gradient at internal- and processor-face
        //  quadrature points
        template<class Type>
        void fGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            CompactListList
            <
                typename outerProduct<vector, Type>::type
            >& result
        ) const;

        //- Disallow default bitwise copy construct
        kExactLeastSquares(const kExactLeastSquares&);

        //- Disallow default bitwise assignment
        void operator=(const kExactLeastSquares&);


public:

    //- Runtime type information
    TypeName("kExactLeastSquares");


    // Constructors

        //- Construct from components
        kExactLeastSquares
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );


    //- Destructor
    virtual ~kExactLeastSquares();


    // Member Functions

        //- Return polynomial order of the local reconstruction
        virtual label polynomialOrder() const
        {
            return polynomialOrder_;
        }

        //- Const access to the weight function
        const weightFunction& weightFunc() const
        {
            if (weightFuncPtr_.empty())
            {
                makeWeightFunction();
            }

            return *weightFuncPtr_;
        }

        //- Return quadrature data
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

        //- Return face-gradient reconstruction coefficients
        virtual const List<CompactListList<vector>>& faceGradCoeffs() const;

        //- Return second-order central cell moments owned by quadrature
        const List<symmTensor>& secondOrderCellMoments() const;

        //- Return third-order central cell moments owned by quadrature
        const List<symmTensor3rdOrder>& thirdOrderCellMoments() const;

        //- Clear owned geometry-dependent data and coefficients
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

        //- Gradient of a scalar field at face quadrature points
        virtual void fGrad
        (
            const volScalarField& vf,
            CompactListList<vector>& result
        ) const;

        //- Gradient of a vector field at face quadrature points
        virtual void fGrad
        (
            const volVectorField& vf,
            CompactListList<tensor>& result
        ) const;

        //- Scalar field values at patch face quadrature points
        virtual autoPtr<CompactListList<scalar>> patchFaceQuadValues
        (
            const volScalarField& vf,
            const label patchI
        ) const;

        //- Vector field values at patch face quadrature points
        virtual autoPtr<CompactListList<vector>> patchFaceQuadValues
        (
            const volVectorField& vf,
            const label patchI
        ) const;

        //- Cell-centre gradient of a scalar field
        virtual tmp<volVectorField> grad(const volScalarField& vf) const;

        //- Cell-centre gradient of a vector field
        virtual tmp<volTensorField> grad(const volVectorField& vf) const;

        //- Cell-centre second gradient of a scalar field
        virtual tmp<volSymmTensorField> secondGrad
        (
            const volScalarField& vf
        ) const;

        //- Cell-centre third gradient of a scalar field
        virtual autoPtr<List<symmTensor3rdOrder>> thirdGrad
        (
            const volScalarField& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kExactLeastSquaresTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
