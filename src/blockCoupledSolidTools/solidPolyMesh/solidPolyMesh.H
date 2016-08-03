/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    solidPolyMesh

Description
    Mesh for solid mechanics allowing implicit cell-point neighbours to enable
    block coupled solution.

    Note on addressing:
    The limitation the fvMesh addressing is that only cell-face neighbours
    are treated implicitly; however, for a fully coupled solution in solid
    mechanics all cell-point neighbours must be implicitly. Therefore, we
    must create a new lduAddressing to allow cell-point neighbours to be
    included in the matrix.

    To achieve this we will define implicit bonds which connect cell
    centres to the cell centres of cell-point neighbours. This way the implicit
    bonds act in a analogous manner to the faces in the fvMesh, where each
    implicit bond has an owner and neighbour. This allows us to create
    lduAddressing in a similar manner to fvMesh; the size of upper and lower
    will be number of implicit bonds. The implicit bonds are created by visiting
    every cellw in ascending order and adding all point neighbour in ascending
    order not already added. This way the starts of the implicit bonds define
    the owners (and lowerAddr) and the ends define the neighbours (and
    upperAddr).

    Notes on boundary faces:
    We now also include the boundary faces implicitly; this is required to
    implement traction boundary conditions implicitly.
    The boundary face unknowns are placed at the end of the list of unknowns;
    potentially this may produce a large bandwidth but it seems to work well.

SourceFiles
    solidPolyMesh.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef solidPolyMesh_H
#define solidPolyMesh_H

#include "GeoMesh.H"
#include "polyMesh.H"
#include "lduMesh.H"
#include "primitiveMesh.H"
#include "solidPolyBoundaryMesh.H"
#include "solidPolyMeshLduAddressing.H"
#include "FieldFields.H"
#include "IOdictionary.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "CompactListList.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidPolyMeshLduAddressing;
class newLeastSquaresVolPointInterpolation;

/*---------------------------------------------------------------------------*\
                     Class solidPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class solidPolyMesh
:
    public GeoMesh<polyMesh>,
    public lduMesh,
    public IOdictionary // so we can look it up
{
    // Permanent data

        //- Store non-const reference to fvMesh
        fvMesh& fMesh_;

        //- Boundary mesh
        solidPolyBoundaryMesh boundary_;


    // Demand-driven data

        //- LDU addressing
        mutable solidPolyMeshLduAddressing* lduPtr_;

        //- Implicit bonds
        mutable edgeList* implicitBondsPtr_;

        //- Whether the implicit bond is through an internal face
        //  These correspond to the neighbours that would have been in the
        //  original fvMesh compact implicit molecule
        mutable boolList* faceImplicitBondsPtr_;

        //- This map tranforms the fvMesh upper/lower to the new addressing
        //  i.e. given an internal face ID, it returns the ID into the upper
        //  and lower
        mutable labelList* fvMeshAddressingMapPtr_;

        //- Cell point cells i.e. cell neighbours connect through any of the
        // cell points
        mutable labelListList* cellPointCellsPtr_;

        //- Cell point boundary faces i.e. boundary face neighbours
        // connected to the cell through any of the face points
        mutable labelListList* cellPointBFacesPtr_;

        //- Boundary face point boundary faces i.e. boundary face neighbours
        // connected to other boundary faces through any of the face points
        mutable labelListList* bFacePointBFacesPtr_;

        //- Cell implicit bonds connected to other cells
        mutable labelListList* cellImplicitBondsPtr_;

        //- Cell implicit bonds connected to boundary faces
        mutable labelListList* cellBFaceImplicitBondsPtr_;

        //- Boundary faces implicit bonds connected to other boundary faces
        mutable labelListList* bFaceImplicitBondsPtr_;

        // Vol-to-point least squares interpolator
        newLeastSquaresVolPointInterpolation volToPointInterp_;

        // Global coefficients
        mutable List< Map<tensorField> >* pointProcFacesCoeffsPtr_;
        mutable List< Map<tensorField> >* pointProcBndFacesCoeffsPtr_;
        mutable List< Map<tensorField> >* pointProcCellsCoeffsPtr_;
        mutable List< Map<tensorField> >* gPtNgbProcBndFaceCoeffsPtr_;
        mutable List< Map<tensorField> >* gPtNgbProcCellCoeffsPtr_;

        // Value of fixed value point components
        mutable Map<vector>* pointFixedComponentPtr_;

        // Direction of fixed value point components
        mutable Map<symmTensor>* pointFixedDirectionPtr_;


    // Private Member Functions

        //- Add global patch. This allows correct matrix multiplications in
        //  parallel
        //void addGlobalPatch();

        //- Disallow default bitwise copy construct
        solidPolyMesh(const solidPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const solidPolyMesh&);


    // Private member functions to calculate demand driven data

        //- Calculate implicit bonds and set faceImplicitBonds
        void calcImplicitBonds() const;

        //- Calculate fvMesh addressing map
        void calcFvMeshAddressingMap() const;

        //- Calculate cell implicit bonds
        void calcCellImplicitBonds() const;

        //- Variable (cell or boundary face) ordering (bandCompression). Adapted
        // Calculate fixed point map fields
        void calcPointFixed(const volVectorField& D) const;

        //- Clear out all demand driven data
        //void clearOut() const;


public:

    // Declare name of the class and its debug switch
    ClassName("solidPolyMesh");

    typedef solidPolyMesh Mesh;
    typedef solidPolyBoundaryMesh BoundaryMesh;

        // Static mesh analysis data

             //- Estimated number of implicit bonds per cell
             static const unsigned implicitBondsPerCell_ = 40;


    // Constructors

        //- Construct from components
        //explicit solidPolyMesh(const polyMesh& pMesh);
        explicit solidPolyMesh(fvMesh&);

    // Destructor

        ~solidPolyMesh();


    // Member Functions

        // Access

            //- Return the top-level database
            const Time& time() const
            {
                return mesh_.time();
            }

            //- Return the object registry
            virtual const objectRegistry& thisDb() const
            {
                return GeoMesh<polyMesh>::thisDb();
            }

            // Return const reference to the volume-to-point interpolator
            const newLeastSquaresVolPointInterpolation& volToPointInterp() const
            {
                return volToPointInterp_;
            }

            //- Return implicit bonds
            const edgeList& implicitBonds() const;

            //- Return near implicit bonds
            const boolList& faceImplicitBonds() const;

            //- Return fvMesh addressing map
            const labelList& fvMeshAddressingMap() const;

            //- Return cell point cells
            const labelListList& cellPointCells() const;

            //- Return cell point boundary faces
            const labelListList& cellPointBFaces() const;

            //- Return boundary face point boundary faces
            const labelListList& bFacePointBFaces() const;

            //- Return cell implicit bonds
            const labelListList& cellImplicitBonds() const;

            //- Return cell implicit bonds connecting with boundary faces
            const labelListList& cellBFaceImplicitBonds() const;

            //- Return boundary face implicit bonds connect to other boundary
            // faces
            const labelListList& bFaceImplicitBonds() const;

            //- Return number of solution variables
            label nVariables() const;

            //- Return reference to boundary mesh
            const solidPolyBoundaryMesh& boundary() const
            {
                 return boundary_;
            }

            //- Return ldu addressing
            const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Owner
            const unallocLabelList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Neighbour
            const unallocLabelList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- Return implicit bond ID connecting two cells
            // Returns -1 if the two cells are not connected
            label findCellCellImplicitBond
            (
                const label firstCellID, const label secondCellID
            ) const;

            //- Return implicit bond ID connecting a cell to a boundary face
            // Returns -1 if they are not connected by an implicit bond
            label findCellFaceImplicitBond
            (
                const label cellID, const label bFaceID
            ) const;

            //- Return implicit bond ID connecting two boundary faces
            // Returns -1 if they are not connected
            label findFaceFaceImplicitBond
            (
                const label firstBFaceID, const label secondBFaceID
            ) const;

            //- Return implicit bond ID connecting the the given IDs
            //  The IDs can be wither the cellID or the boundaryFaceID
            label findBond
            (
                const label firstID, const label secondID
            ) const;

            //- Return the variable ID (index in variable list) given the faceID
            //  Returns -1 for internal faces and boundary faces which are empty
            //  or processor boundaries
            //  Note: initially the varibales are ordered as the cells followed
            //  by boundary faces: findOldVariableID returns the index into this
            //  initial ordering. This list is then renumbered (band
            //  compression): the findNewVariableID function returns the index
            //  into the renumbered list
            label findOldVariableID(const label bFaceID) const;

            //- Returns the variable ID given the cellID or boundary face ID
            //  Returns -1 for invalid IDs
            //  See note in above findOldVariableID
            label findNewVariableID(const label cellOrBFaceID) const;

            //- Return map from old variable IDs to new variable IDs
            const labelList& oldToNewAddressing() const;

            //- Add the given fvVectorMatrix to the given block matrix and block
            //  source
            void addFvMatrix
            (
                BlockLduMatrix<vector>& blockM,
                vectorField& blockB,
                const fvVectorMatrix& fvM,
                const bool flipSign = false
            ) const;

            //- Add the given fvVectorMatrix source to the given block source
            void addFvSource
            (
                vectorField& blockB,
                const volVectorField& fvB
            ) const;

            //- Insert the boundary condition equations
            void insertBoundaryConditions
            (
                BlockLduMatrix<vector>& blockM,
                vectorField& blockB,
                const surfaceScalarField& muf,
                const surfaceScalarField& lambdaf,
                const volVectorField& D
            ) const;


            // Parallel mesh analysis data

                //- Return parallel info
                const globalMeshData& globalData() const
                {
                    return mesh_.globalData();
                }

                //- Make global coefficients
                void makeGlobalCoeffs
                (
                    const newLeastSquaresVolPointInterpolation& volToPointInterp
                ) const;

                //- Global coefficients
                List< Map<tensorField> >& pointProcFacesCoeffs() const;
                List< Map<tensorField> >& pointProcBndFacesCoeffs() const;
                List< Map<tensorField> >& pointProcCellsCoeffs() const;
                List< Map<tensorField> >& gPtNgbProcBndFaceCoeffs() const;
                List< Map<tensorField> >& gPtNgbProcCellCoeffs() const;

                //- Update global fields
                // This is called from Amul and calculates the multiplication
                // for the pointCells of global shared points
                void updateGlobalFields
                (
                    vectorField& Ax,
                    const vectorField& x
                ) const;

                //- Copy the given solution vector into the volVectorField
                void copySolutionVector
                (
                    const vectorField& x,
                    volVectorField& vf
                ) const;

            // Return fixed point components map
            const Map<vector>& pointFixedComponent
            (
                const volVectorField& D
            ) const;

            // Return fixed point directions map
            const Map<symmTensor>& pointFixedDirection
            (
                const volVectorField& D
            ) const;

        // Edit

            //- Clear out all demand driven data
            void clearOut() const;

            //- Clear out global coeffs
            void clearOutGlobalCoeffs() const;


    // Member Operators

        bool operator!=(const solidPolyMesh&) const;
        bool operator==(const solidPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
