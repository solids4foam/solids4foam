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
    Foam::dualMechanicalModel

Description
    This is a stripped down version of the mechancialModel class, which is
    used by the vertex-centred finite volume solid mechanics solvers.

    This class uses the mechancialProperties dict from a mechancialModel object
    (e.g. associated with the primary mesh) to create its own mechanical laws
    that are associated with another mesh (e.g. the dual mesh).

SourceFiles
    dualMechanicalModel.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef dualMechanicalModel_H
#define dualMechanicalModel_H

#include "mechanicalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class dualMechanicalModel Declaration
\*---------------------------------------------------------------------------*/

class dualMechanicalModel
:
    public PtrList<mechanicalLaw>
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Const reference to the primary mesh mechanical model
        const mechanicalModel& mechModel_;

        //- Is the current case plane stress or plain strain/3-D
        //  Todo: this should be moved to actual mechanical law
        const Switch planeStress_;

        //- Does the solidModel use an incremental approach?
        //  i.e. does the solidModel solve for DD instead of D
        const bool incremental_;

        //- Copy of dualFaceToCellMap list
        //  We take a copy to avoid an possible memory issues if the field is
        //  freed
        const labelList dualFaceToCell_;

        //- For each material, return a field of 1s and 0s, where primary mesh
        //  cells in this material are set to 1 and the rest of the faces are
        //  set to 0
        mutable PtrList<volScalarField> cellInThisMaterialList_;

        //- For each material, return a field of 1s and 0s, where dual faces in
        //  this material are set to 1 and the rest of the faces are set to 0
        mutable PtrList<surfaceScalarField> dualFaceInThisMaterialList_;

    // Private Member Functions

        //- Make cellInThisMaterialList
        void makeCellInThisMaterialList() const;

        //- Return const acess to cellInThisMaterialList
        const PtrList<volScalarField>& cellInThisMaterialList() const;

        //- Make dualFaceInThisMaterialList
        void makeDualFaceInThisMaterialList() const;

        //- Return const acess to dualFaceInThisMaterialList
        const PtrList<surfaceScalarField>& dualFaceInThisMaterialList() const;

        //- Clear out demand driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        dualMechanicalModel(const dualMechanicalModel&);

        //- Disallow default bitwise assignment
        void operator=(const dualMechanicalModel&);


public:

    // Constructors

        //- Construct from an fvMesh
        dualMechanicalModel
        (
            const fvMesh& mesh,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const bool incremental,
            const mechanicalModel& mechModel,
            const labelList& dualFaceToCell
        );


    //- Destructor
    ~dualMechanicalModel();


    // Member Functions

        // Access

            //- Return const access to the mesh
            const fvMesh& mesh() const;

#ifdef OPENFOAM_NOT_EXTEND
            //- Material tangent field defined at dual mesh faces
            tmp< Field< RectangularMatrix<scalar>>>
            materialTangentFaceField() const;
#endif

            //- Implicit modulus field at the faces
            tmp<surfaceScalarField> impKf() const;

            //- Bulk modulus field at the faces
            tmp<surfaceScalarField> bulkModulus() const;

        // Edit

            //- Calculate the stress surfaceField
            void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
