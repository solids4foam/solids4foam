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
    higherOrderGrad

Description
    Class for calculating higher order cell-centre gradients using the least
    squares method.

    Currently, linear and quadratic polynomials are implemented.

    Cells are added to the stencil in layers.

Author
    Ivan Batistic
    Philip Cardiff

SourceFiles
    higherOrderGrad.C

\*---------------------------------------------------------------------------*/

#ifndef higherOrderGrad_H
#define higherOrderGrad_H

#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class higherOrderGrad Declaration
\*---------------------------------------------------------------------------*/

class higherOrderGrad
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Order of interpolation
        const label N_;

        //- Number of layers in the stencils
        const label nLayers_;

        //- Kernel shape parameter
        const label k_;

        //- Maxium stencil size: used for initialising dynamic lists
        const label maxStencilSize_;

        //- Calculate the condition number field
        const bool calcConditionNumber_;

        //- Condition number field
        mutable autoPtr<volScalarField> condNumberPtr_;

        //- Least squares cell stencils
        mutable autoPtr<List<DynamicList<label>>> stencilsPtr_;

        //- Interpolation coefficients
        mutable autoPtr<List<DynamicList<scalar>>> interpCoeffsPtr_;

        //- Gradient coefficients
        mutable autoPtr<List<DynamicList<vector>>> interpGradCoeffsPtr_;

        //- Cell boundary face indices
        mutable autoPtr<List<DynamicList<label>>> cellBoundaryFacesPtr_;


    // Private Member Functions

        //- Least squares cell stencils
        void makeStencils() const;

        //- Return the least squares cell stencils
        const List<DynamicList<label>> stencils() const;

        //- Calculate interpolation and gradient coefficients
        void calcCoeffs() const;

        //- Return the interpolation coefficients
        const List<DynamicList<scalar>>& interpCoeffs() const;

        //- Return the gradient coefficients
        const List<DynamicList<vector>>& interpGradCoeffs() const;

        //- Return the cell boundary face indices
        const List<DynamicList<label>>& cellBoundaryFaces() const;

        //- Disallow default bitwise copy construct
        higherOrderGrad(const higherOrderGrad&);

        //- Disallow default bitwise assignment
        void operator=(const higherOrderGrad&);


public:

    //- Runtime type information
    TypeName("higherOrderGrad");

    // Constructors

        //- Construct from components
        higherOrderGrad
        (
            const fvMesh& mesh,
            const label N = 2,
            const label nLayers = 4,
            const label k = 6,
            const label maxStencilSize = 70,
            const bool calcConditionNumber = false
        );

    // Destructor

        virtual ~higherOrderGrad();


    // Member Functions

        // Access

            //- Condition number field
            const volScalarField& condNumber() const
            {
                if (calcConditionNumber_)
                {
                    FatalErrorInFunction
                        << "Field not calculated!" << abort(FatalError);
                }

                return *condNumberPtr_;
            }

        // Edit

            //- Calculate the gradient of a volVectorField
            tmp<volTensorField> grad(const volVectorField& D);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
