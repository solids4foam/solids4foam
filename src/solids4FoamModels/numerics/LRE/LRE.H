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
    LRE

Description
    The LRE method is used to estimate the derivatives of an interpolation
    function. It provides a truncated Taylor expansion of the function
    by approximating the derivatives. The Taylor series includes terms for
    the function value and its partial derivatives with respect
    to each variable.

    Cells are added to the stencil in layers. Face stencil is build using
    face cells stencils.

    Source:
        Castrillo P, Canelas A, Schillaci E, Rigola J, Oliva A. High-order
	finite volume method for linear elasticity on unstructured meshes.
	Comput Struct, 2022.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    LRE.C

\*---------------------------------------------------------------------------*/

#ifndef LRE_H
#define LRE_H

#include "fvMesh.H"
#include "globalIndex.H"
#include <Eigen/Dense>
#include "tensorList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class LRE Declaration
\*---------------------------------------------------------------------------*/

class LRE
{
public:
        //- Weight functions
        enum weightFunction
	{
	    ONE,
	    RAD_SYMM_EXP
	};

	//- Weight funtions names
        static const Enum<weightFunction> weightFunctionNames_;

private:
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Include patch in the stencils
        const boolList includePatchInStencils_;

        //- Order of interpolation
        const label N_;

        //- Number of layers in the stencils
        const label nLayers_;

        //- Kernel shape parameter
        const scalar k_;

        //- Weight function
        const weightFunction weightFunction_;

        //- Maxium stencil size: used for initialising dynamic lists
        const label maxStencilSize_;

        //- Global cell addressing object for converting between local and
        //  global indices
        const globalIndex globalCells_;

        //- Use QR decomposition to solve the least squares problem; if false,
        //  use Cholesky decomposition to solve the normal equations
        //  QR is more stable (to round-off error) but more expensive
        const Switch useQRDecomposition_;

        //- Use new parallelised global cell stencils approach
        //  Once this has been verified, we can remove this flag and make this
        //  the default method
        const Switch useGlobalStencils_;

        //- Calculate the condition number field
        const Switch calcConditionNumber_;

        //- Condition number field of the least squares system
        mutable autoPtr<volScalarField> conditionNumberPtr_;

        //- Least squares cell stencils
        mutable autoPtr<List<DynamicList<label>>> stencilsPtr_;

        //- Cell boundary face indices
        mutable autoPtr<List<DynamicList<label>>> stencilsBoundaryFacesPtr_;

        //- Least squares global cell stencils
        //  This is a labelList for each cell on this processor, containing the
        //  global cell indices of all cells in its stencil, including those on
        //  other processors
        //  The globalCells object can be used to check if a cell is local or
        //  not and to convert between global and local indices
        mutable autoPtr<labelListList> globalCellStencilsPtr_;

        //- Least squares global face stencils
        //  This is a labelList for each face (internal and boundary) on this
        //  processor, containing the global cell indices of all cells in its
        //  stencil, including those on other processors
        //  The globalCells object can be used to check if a cell is local or
        //  not and to convert between global and local indices
        mutable autoPtr<labelListList> globalFaceStencilsPtr_;

        //- QR Interpolation coefficients
        mutable autoPtr<List<DynamicList<scalar>>> QRInterpCoeffsPtr_;

        //- QR Gradient cell coefficients
        mutable autoPtr<List<DynamicList<vector>>> QRGradCoeffsPtr_;

        //- QR Gradient face coefficients
        mutable autoPtr<List<DynamicList<vector>>> QRGradFaceCoeffsPtr_;

        //- QR Gradient face Gauss points coefficients
        mutable autoPtr<List<List<DynamicList<vector>>>> QRGradFaceGPCoeffsPtr_;

        //- Cholesky decomposition
        mutable autoPtr<List<Eigen::LLT<Eigen::MatrixXd>>> choleskyPtr_;

        //- Qhat matrix needed for Cholesky decomposition
        mutable autoPtr<List<Eigen::MatrixXd>> QhatPtr_;

        //- sqrtW matrix needed for Cholesky decomposition
        mutable autoPtr
        <
            List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>
        > sqrtWPtr_;

        //- Face Gauss points
        mutable autoPtr<List<List<point>>> faceQuadPointsPtr_;

        //- Face Gauss points weights
        mutable autoPtr<List<List<scalar>>> faceQuadPointsWeightPtr_;

    // Private Member Functions

        //- Make the least squares cell stencils
        void makeGlobalCellStencils() const;

        //- Make the least squares face stencils
        void makeGlobalFaceStencils() const;

        //- Least squares cell stencils
        void makeStencils() const;

        //- Return the least squares cell stencils
        const List<DynamicList<label>>& stencils() const;

        //- Return the least squares global cell stencils
        const labelListList& globalCellStencils() const;

        //- Weight using radially symmetric exponential function
        scalar radiallySymmetricExponentialWeight
	(
	     const scalar d,
	     const scalar maxDist
	) const;

        //- Weight funtcion
        scalar weight(const scalar d, const scalar maxDist) const;

        //- Calculate exponents for the Tayler series expansion terms
        void generateExponents
        (
            const label N,
            DynamicList<FixedList<label, 3>>& exponents
        ) const;

        //- Calculate coefficients for the QR decomposition approach
        void calcQRCoeffs() const;

        //- Calculate coefficients for the QR decomposition approach using the
        //  the global cell stencils
        void calcGlobalQRCoeffs() const;

        //- Calculate coefficients for the QR decomposition approach using the
        //  global cell stencils for the faces
        void calcGlobalQRFaceCoeffs() const;

        //- Calculate coefficients for the QR decomposition approach at face
        //  Gauss points using the global cell stencils for the faces
        void calcGlobalQRFaceGPCoeffs() const;

        //- Calculate coefficients for the Cholesky decomposition approach
        void calcCholeskyCoeffs() const;

        //- Calculate coefficients for the Cholesky decomposition approach using
        //  the global cell stencils
        void calcGlobalCholeskyCoeffs() const;

        //- Calculate Gauss points and corresponding weights in 3D
        void calcQuadPointsAndWeights3D() const;

        //- Calculate Gauss points and corresponding weights in 2D
        void calcQuadPointsAndWeights2D() const;

        //- Calls the same funtion implemented for 3D or for 2D meshes
        void calcQuadPointsAndWeights() const;

        //- Return the interpolation coefficients
        const List<DynamicList<scalar>>& QRInterpCoeffs() const;

        //- Return the gradient coefficients
        const List<DynamicList<vector>>& QRGradCoeffs() const;

        //- Return the gradient face coefficients
        const List<DynamicList<vector>>& QRGradFaceCoeffs() const;

        //- Return the cell boundary face indices
        const List<DynamicList<label>>& stencilsBoundaryFaces() const;

        //- Return Cholesky matrices
        const List<Eigen::LLT<Eigen::MatrixXd>>& cholesky() const;

        //- Return Qhat matirces used by the Cholesky decomposition approach
        const List<Eigen::MatrixXd>& Qhat() const;

        //- Return sqrtW weight matirces used by the Cholesky decomposition
        //  approach
        const List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>&
            sqrtW() const;

        //- Return the condition number field
        volScalarField& conditionNumber() const;

        //- Make the condition number field
        void makeConditionNumber() const;

        //- Request data for off-processor cells in the stencils
        //  Given the local processor field data, return a map containing
        //  the field data from other processors that is required for the local
        //  processor cell stencils
        template<class Type>
        void requestGlobalStencilData
        (
            const Field<Type>& localField,
            Map<Type>& globalField
        ) const;


        //- Disallow default bitwise copy construct
        LRE(const LRE&);

        //- Disallow default bitwise assignment
        void operator=(const LRE&);


public:

    //- Runtime type information
    TypeName("LRE");

    // Constructors

        //- Construct from components
        LRE
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );

    // Destructor

        virtual ~LRE();


    // Member Functions

            //- Const access tot he weight function enum
            const weightFunction& weightFunc() const
            {
		return weightFunction_;
	    }

            //- Return order of accuracy
            label order() const {return N_;};

            //- Return overall cell stencil size using cell faces stencils
            tmp<labelField> cellFacesStencilSize() const;

            //- Return the face Gauss points weights
            const List<List<scalar>>& faceGaussPointsWeight() const;

            //- Return the face Gauss points
            const List<List<point>>& faceGaussPoints() const;

            //- Return the least squares face cell stencils
            const labelListList& globalFaceStencils() const;

            //- Return the gradient face Gauss points coefficients
            const List<List<DynamicList<vector>>>& QRGradFaceGPCoeffs() const;

        // Edit

            //- Calculate the gradient of a volVectorField
            tmp<volTensorField> grad(const volVectorField& D) const;

            //- Calculate the face gradient of a volVectorField
            tmp<surfaceTensorField> fGrad(const volVectorField& D) const;

            //- Calculate the face Gauss points gradient of a volVectorField
            autoPtr<List<List<tensor>>> gradDQuad
	    (
	        const volVectorField& D
	    ) const;

            //- Calculate the gradient of a volVectorField using QR
            //  decomposition
            tmp<volTensorField> gradQR(const volVectorField& D) const;

            //- Calculate the gradient of a volVectorField using QR
            //  decomposition using the global cell stencils
            tmp<volTensorField> gradGlobalQR(const volVectorField& D) const;

            //- Calculate the gradient of a volVectorField using Cholesky
            //  decomposition
            tmp<volTensorField> gradCholesky(const volVectorField& D) const;

            //- Calculate the gradient of a volVectorField using Cholesky
            //  decomposition using the global cell stencils
            tmp<volTensorField> gradGlobalCholesky(const volVectorField& D) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LRETemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
