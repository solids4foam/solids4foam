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


private:

    // Private Member Functions

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

        //- Clear owned geometry-dependent data and coefficients
        virtual void clear() const = 0;


    // Evaluation

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
