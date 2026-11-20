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
    foamPetscMeshAssembly

Description
    Assembly of finite volume discretisations on one mesh into PETSc matrices
    and vectors.

    The assembly holds the local-to-global numbering of the unknowns of the
    mesh, which are located at the cells or at the points, and the data of the
    neighbouring processor cells. It provides the preallocation of a Jacobian
    matrix and a solution vector for a given number of unknowns per cell or
    point, the insertion of fvMatrix and other operator coefficients into a
    Jacobian, and the transfer of fields to and from PETSc vectors.

    One assembly is used per mesh. A monolithic system of several meshes uses
    the assembly of each mesh for its diagonal block.

Author
    Philip Cardiff, UCD.

SourceFiles
    foamPetscMeshAssembly.C
    foamPetscMeshAssemblyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamPetscMeshAssembly_H
#define foamPetscMeshAssembly_H

#ifdef USE_PETSC

#include <petscsnes.h>
#include "globalIndex.H"
#include "globalPointIndices.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "leastSquaresS4fVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class foamPetscMeshAssembly Declaration
\*---------------------------------------------------------------------------*/

class foamPetscMeshAssembly
{
public:

        //- Location of the solution unknowns
        enum solutionLocation
        {
            CELLS,
            POINTS,
            NONE
        };

        //- Solution location names
#ifdef OPENFOAM_COM
        static const Enum<solutionLocation> solutionLocationNames_;
#else
        static const NamedEnum<solutionLocation, 3> solutionLocationNames_;
#endif


private:

    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Location of the solution unknowns, e.g. cells or points
        const solutionLocation location_;

        //- Helper object for converting between local (on processor) and global
        //  cell indices
        autoPtr<globalIndex> globalCellsPtr_;

        //- Helper object for converting between local (on processor) and global
        //  point indices
        autoPtr<globalPointIndices> globalPointsPtr_;

        //- Global cell IDs for the neighbouring processor patch cells
        mutable PtrList<labelList> neiProcGlobalIDs_;

        //- Cell volumes for the neighbouring processor patch cells
        mutable PtrList<scalarField> neiProcVolumes_;


    // Private Member Functions

        //- Check that the given mesh is the mesh of this assembly
        void checkMesh(const fvMesh& mesh) const;

        //- Make the neighbour processor patch fields for the given mesh
        void makeNeiProcFields(const fvMesh& mesh) const;

        //- Return a reference to neiProcGlobalIDs for the given mesh
        const PtrList<labelList>& neiProcGlobalIDs(const fvMesh& mesh) const;

        //- Return a reference to neiProcGlobalIDs for the given mesh
        const PtrList<scalarField>& neiProcVolumes(const fvMesh& mesh) const;

        //- Return a reference to the least squares vectors
        const leastSquaresS4fVectors& lsVectors(const volScalarField& p) const;

        //- Disallow default bitwise copy construct
        foamPetscMeshAssembly(const foamPetscMeshAssembly&);

        //- Disallow default bitwise assignment
        void operator=(const foamPetscMeshAssembly&);


public:

    //- Runtime type information
    TypeName("foamPetscMeshAssembly");


    // Constructors

        //- Construct from the mesh and the location of the unknowns
        foamPetscMeshAssembly
        (
            const fvMesh& mesh,
            const solutionLocation location = solutionLocation::CELLS
        );


    //- Destructor
    virtual ~foamPetscMeshAssembly()
    {}


    // Member Functions

        // Access

            //- Return the mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return the location of the unknowns
            solutionLocation location() const
            {
                return location_;
            }

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                if (globalCellsPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not set" << exit(FatalError);
                }

                return globalCellsPtr_();
            }

            //- Const access to the global point object
            const globalPointIndices& globalPoints() const
            {
                if (globalPointsPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not set" << exit(FatalError);
                }

                return globalPointsPtr_();
            }


        // Edit

            //- Initialise a Jacobian with the compact stencil of the mesh and
            //  the given blockSize (number of scalar equations per cell or
            //  point)
            label initialiseJacobian
            (
                Mat& jac,
                const label blockSize,
                const bool createMat = true  // MatCreate or similar should be called
            ) const;

            //- Initialise a solution vector for the mesh and the given
            //  blockSize
            label initialiseSolution
            (
                Vec& x,
                const label blockSize,
                const bool createVec = true  // VecCreate or similar should be called
            ) const;

            //- Inserts the OpenFOAM fvMatrix 'fvM' into the PETSc matrix 'jac'.
            //  Converts the finite volume matrix into block coefficients using
            //  the specified row and column offsets (rowOffset, colOffset)
            //  within each block.
            //  nScalarEqns indicates the number of scalar equations to add, e.g.
            //  3 for a 3-D vector equation, 2 for 2-D vector equation, 1 for a
            //  scalar equation
            //  A non-zero return value indicates an error.
            //  Note: fvM should be const but OF9 does not have required const
            //  access functions
            template <class Type>
            label InsertFvMatrixIntoPETScMatrix
            (
                fvMatrix<Type>& fvM,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns
            ) const;

            //- Inserts least squares gradient into a PETSc matrix at location
            //  (rowOffset, colOffset) within the block
            //  Only implemented for the gradient of a scalar, e.g. grad(p)
            //  If "flipSign" is true, -grad(p) is inserted rather than grad(p)
            label InsertFvmGradIntoPETScMatrix
            (
                const volScalarField& p,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = -1.0
            ) const;

            //- Insert V*div(gamma*SfHat & interpolate(grad(p))) into a PETSc
            //  matrix. This is the least-squares gradient-interpolation part of
            //  the diff-stencil pressure stabilisation.
            label InsertFvcDivGradInterpolateIntoPETScMatrix
            (
                const volScalarField& p,
                const surfaceScalarField& gamma,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const scalar scale = 1.0
            ) const;

            //- Insert the missing tensor part of the linearisation of
            //  div(phi,U) into a PETSc matrix, where phi = Sf & Uf.
            //  The caller must insert the segregated fvm::div(phi,U)
            //  contribution separately. An upwind discretisation is assumed.
            label InsertFvmDivPhiUIntoPETScMatrix
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = -1.0
            ) const;

            //- Insert div(U) into the scalar p equation of a PETSc matrix.
            //  The inserted matrix is `scale * V * div(U)`. Callers
            //  needing `-V * div(U)` should pass `scale = -1.0`, or
            //  `-alpha` if they additionally want to row-scale by `alpha`
            //  (e.g. for a non-dimensionalised pressure residual).
            //  Includes proper processor-patch coupling and a Dirichlet
            //  boundary correction that undoes the implicit zeroGradient
            //  assumption baked into negSumDiag at boundary cells.
            label InsertFvmDivUIntoPETScMatrix
            (
                const volScalarField& p,
                const volVectorField& U,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = 1.0
            ) const;

            //- Insert the Gauss-style grad(p) into the vector D equation
            //  of a PETSc matrix. The inserted matrix is
            //  `scale * V * grad(p)`. Callers needing `-V*grad(p)`
            //  should pass `scale = -1.0`. This is the structural
            //  transpose of `InsertFvmDivUIntoPETScMatrix` (same per-face
            //  coefficient pattern, inserted at block position
            //  (rowOffset + cmptI, colOffset) instead of
            //  (rowOffset, colOffset + cmptI)).
            label InsertFvmGradPGaussIntoPETScMatrix
            (
                const volScalarField& p,
                Mat jac,
                const label rowOffset,
                const label colOffset,
                const label nScalarEqns,
                const scalar scale = 1.0
            ) const;

            //- Extract a subset of components a raw scalar array 'x' into the
            //  OpenFOAM field 'vf'. The array 'x' is organized in blocks of
            //  'xBlockSize' scalars per block, and extraction begins at the
            //  specified 'offset'. For each element in 'vf', the function
            //  inserts  the components specified in 'compIndices'. If
            //  'compIndices' is empty, it copies all available components
            //  (i.e., the first pTraits<Type>::nComponents scalars), assuming
            //  that 'vf' is stored contiguously.
            //
            //  Parameters:
            //    vf          - The OpenFOAM field (e.g., Field<scalar>,
            //                  Field<vector>, Field<tensor>) into which the
            //                  components will be inserted.
            //    x           - Pointer to the raw scalar array from where the
            //                  components will be extracted.
            //    offset      - Starting index within each block in 'x' where
            //                  extraction begins.
            //    xBlockSize  - The total number of scalar values stored per
            //                  block in 'x' (i.e., the stride).
            //    compIndices - (Optional) List of component indices specifying
            //                  which scalars from each element of 'vf' should
            //                  be inserted. If empty, all components are copied.
            template <class Type>
            void ExtractFieldComponents
            (
                const PetscScalar *x,
                Field<Type>& vf,
                const label offset,
                const label xBlockSize,
                const labelList& compIndices = labelList()
            ) const;

            //- Same as above, but takes a PETSc "Vec x" instead of a raw
            //  PetscScalar array pointer. The xBlockSize is looked up
            //  from "Vec x"
            //  Do NOT use this function if x has a non-constant block size
            //  Use the version of this function above, where the block size
            //  is explicitly passed
            template <class Type>
            void ExtractFieldComponents
            (
                const Vec x,
                Field<Type>& vf,
                const label offset,
                const labelList& compIndices = labelList()
            ) const;

            //- Inserts a subset of components from the OpenFOAM field 'vf' into
            //  a raw scalar array 'x'. The array 'x' is organized in blocks of
            //  'blockSize' scalars per cell, and insertion begins at the
            //  specified 'offset'. For each element in 'vf', the function
            //  copies the components specified in 'compIndices'. If
            //  'compIndices' is empty, it copies all available components
            //  (i.e., the first pTraits<Type>::nComponents scalars), assuming
            //  that 'vf' is stored contiguously.
            //
            //  Parameters:
            //    vf          - The OpenFOAM field (e.g., Field<scalar>,
            //                  Field<vector>, Field<tensor>) from which the
            //                  components are extracted.
            //    x           - Pointer to the raw scalar array where the
            //                  components will be inserted.
            //    offset      - Starting index within each block in 'x' where
            //                  insertion begins.
            //    xBlockSize  - The total number of scalar values stored per
            //                  block in 'x' (i.e., the stride).
            //    compIndices - (Optional) List of component indices specifying
            //                  which scalars from each element of 'vf' should
            //                  be copied. If empty, all components are copied.
            template <class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                PetscScalar *x,
                const label offset,
                const label xBlockSize,
                const labelList& compIndices = labelList()
            ) const;

            //- Same as above, but takes a PETSc "Vec x" instead of a raw
            //  PetscScalar array pointer. The xBlockSize is looked up
            //  from "Vec x"
            //  Do NOT use this function if x has a non-constant block size
            //  Use the version of this function above, where the block size
            //  is explicitly passed
            template <class Type>
            void InsertFieldComponents
            (
                const Field<Type>& vf,
                Vec x,
                const label offset,
                const labelList& compIndices = labelList()
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "foamPetscMeshAssemblyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // #ifdef USE_PETSC

#endif

// ************************************************************************* //
