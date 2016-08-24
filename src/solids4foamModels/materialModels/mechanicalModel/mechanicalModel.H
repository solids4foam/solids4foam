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
    Foam::mechanicalModel

Description
    The mechanical model class takes care of reading the mechanical properties
    and creating mechanical laws.

    If there are multiple mechanical laws then cellZones are read, and
    corrections are applied at bi-material interfaces to ensure continuity of
    stress.

    Multiple mechanical laws are treated by defining a separate sub-mesh for
    each mechanical law region (cellZone).

SourceFiles
    mechanicalModel.C
    mechanicalModelTemplates.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalModel_H
#define mechanicalModel_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "mechanicalLaw.H"
#include "Switch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "newFvMeshSubset.H"
#include "newLeastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class mechanicalModel Declaration
\*---------------------------------------------------------------------------*/

class mechanicalModel
:
    public IOdictionary,
    public PtrList<mechanicalLaw>
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Is the current case plane stress or plain strain/3-D
        const Switch planeStress_;

        //- Name of the cellZone for each material
        mutable wordList cellZoneNames_;

        //- Sub-mesh for each material
        //  Not used if there is only one material
        mutable PtrList<newFvMeshSubset> subMeshes_;

        //- Volume field to point field interpolator for the base mesh
        mutable newLeastSquaresVolPointInterpolation* volToPointPtr_;

        //- Volume field to point field interpolator for each mechanical law
        mutable PtrList<newLeastSquaresVolPointInterpolation>
            subMeshVolToPoint_;

        //- Stress volField for each mechanical law
        mutable PtrList<volSymmTensorField> subMeshSigma_;

        //- Stress surfaceField for each mechanical law
        mutable PtrList<surfaceSymmTensorField> subMeshSigmaf_;

        //- Displacement volField for each mechanical law
        mutable PtrList<volVectorField> subMeshD_;

        //- Displacement gradient volField for each mechanical law
        mutable PtrList<volTensorField> subMeshGradD_;

        //- Displacement gradient surfaceField for each mechanical law
        mutable PtrList<surfaceTensorField> subMeshGradDf_;

        //- Displacement pointField for each mechanical law
        mutable PtrList<pointVectorField> subMeshPointD_;

        //- List of faces in the base mesh that are on a bi-material interface
        //  ie. they are shared by two materials
        mutable labelList* interfaceBaseFacesPtr_;

        //- Index of the subMesh on the other side of a bi-material interface
        mutable PtrList<labelList> interfaceShadowSubMeshID_;

        //- Index of the subMesh on the other side of a bi-material interface
        mutable PtrList<labelList> interfaceShadowPatchID_;

        //- Index of the subMesh on the other side of a bi-material interface
        mutable PtrList<labelList> interfaceShadowFaceID_;

        //- The implicit stiffness surface field for Rhie-Chow correction
        mutable surfaceScalarField* impKfcorrPtr_;

        //- Number of materials adjacent to each point in the base mesh
        mutable labelList* pointNumOfMaterialsPtr_;

        //- Isolated interface points
        mutable labelList* isolatedInterfacePointsPtr_;


    // Private Member Functions

        //- Make material sub-meshes
        void makeSubMeshes() const;

        //- Return a const reference to the subMeshes
        const PtrList<newFvMeshSubset>& subMeshes() const;

        //- Make base mesh vol-to-point interpolator
        void makeVolToPoint() const;

        //- Make subMesh vol-to-point interpolators
        void makeSubMeshVolToPoint() const;

        //- Return a const reference to the vol-to-point interpolators
        const PtrList<newLeastSquaresVolPointInterpolation>&
            subMeshVolToPoint() const;

        //- Check the cell zones are defined for each material
        void checkCellZones() const;

        //- Calculate the subMesh sigma fields
        void calcSubMeshSigma() const;

        //- Return the subMesh sigma volFields
        PtrList<volSymmTensorField>& subMeshSigma();

        //- Return the subMesh sigma volFields
        const PtrList<volSymmTensorField>& subMeshSigma() const;

        //- Calculate the subMesh sigmaf fields
        void calcSubMeshSigmaf() const;

        //- Return the subMesh sigma surfaceFields
        PtrList<surfaceSymmTensorField>& subMeshSigmaf();

        //- Return the subMesh sigma surfaceFields
        const PtrList<surfaceSymmTensorField>& subMeshSigmaf() const;

        //- Calculate the subMesh D fields
        void calcSubMeshD() const;

        //- Return the subMesh D volFields
        PtrList<volVectorField>& subMeshD();

        //- Return the subMesh D volFields
        const PtrList<volVectorField>& subMeshD() const;

        //- Calculate the subMesh gradD fields
        void calcSubMeshGradD() const;

        //- Return the subMesh gradD volFields
        PtrList<volTensorField>& subMeshGradD();

        //- Return the subMesh gradD volFields
        const PtrList<volTensorField>& subMeshGradD() const;

        //- Calculate the subMesh gradDf fields
        void calcSubMeshGradDf() const;

        //- Return the subMesh gradDf volFields
        PtrList<surfaceTensorField>& subMeshGradDf();

        //- Return the subMesh gradDf volFields
        const PtrList<surfaceTensorField>& subMeshGradDf() const;

        //- Calculate the subMesh pointD fields
        void calcSubMeshPointD() const;

        //- Return the subMesh pointD volFields
        PtrList<pointVectorField>& subMeshPointD();

        //- Return the subMesh pointD volFields
        const PtrList<pointVectorField>& subMeshPointD() const;

        //- Make the interface base faces
        void makeInterfaceBaseFaces() const;

        //- Return the interface base faces
        const labelList& interfaceBaseFaces() const;

        //- Calculate the interface shadow subMesh indices
        void calcInterfaceShadowIDs() const;

        //- Return the interface shadow subMesh indices
        const PtrList<labelList>& interfaceShadowSubMeshID() const;

        //- Return the interface shadow patch indices
        const PtrList<labelList>& interfaceShadowPatchID() const;

        //- Return the interface shadow face indices
        const PtrList<labelList>& interfaceShadowFaceID() const;

        //- Calculate the impKf field for Rhie-Chow correction
        //  This field is zero on bi-material interfaces
        void calcImpKfcorr() const;

        //- Make materials number for points
        void makePointNumOfMaterials() const;

        //- Return material number for points
        const labelList& pointNumOfMaterials() const;

        //- Make isolated interface points
        void makeIsolatedInterfacePoints() const;

        //- Return isolated interface points
        const labelList& isolatedInterfacePoints() const;

        //- Interpolate the base D to the subMesh D, where we apply corrections
        //  on bi-material interfaces
        void interpolateDtoSubMeshD
        (
            const volVectorField& D,
            const bool useVolFieldSigma
        );

        //- Correct snGrad on bi-material interface patch for subMeshGradD field
        void correctInterfaceSnGrad();

        //- Correct snGrad on bi-material interface patch for subMeshGradDf
        //  field
        void correctInterfaceSnGradf();

        //- Return the implicit stiffness surface field for Rhie-Chow correction
        const surfaceScalarField& impKfcorr() const;

        //- Map a volField from the subMesh to the base mesh
        template<class Type>
        void mapSubMeshVolFields
        (
            const PtrList<GeometricField<Type, fvPatchField, volMesh> >&
                subMeshFields,
            GeometricField<Type, fvPatchField, volMesh>& baseMeshField
        ) const;

        //- Map a surfaceField from the subMesh to the base mesh
        template<class Type>
        void mapSubMeshSurfaceFields
        (
            const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >&
                subMeshFields,
            GeometricField<Type, fvsPatchField, surfaceMesh>& baseMeshField
        ) const;

        //- Map pointFields from the subMesh to the base mesh
        template<class Type>
        void mapSubMeshPointFields
        (
            const PtrList<GeometricField<Type, pointPatchField, pointMesh> >&
                subMshFields,
            GeometricField<Type, pointPatchField, pointMesh>& baseMeshField
        ) const;

        //- Clear out demand driven data
        void clearOut() const;

        //- Disallow default bitwise copy construct
        mechanicalModel(const mechanicalModel&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalModel&);


public:

    // Constructors

        //- Construct from an fvMesh
        mechanicalModel(const fvMesh& mesh);


    //- Destructor
    ~mechanicalModel();


    // Member Functions

        // Access

            //- Return const access to the mesh
            const fvMesh& mesh() const;

            //- Return a const reference to the vol-to-point base mesh
            //  interpolator
            const newLeastSquaresVolPointInterpolation& volToPoint() const;

            //- Return initial density
            tmp<volScalarField> rho() const;

            //- Return the implicit stiffness
            //  This is the diffusivity for the Laplacian term
            tmp<volScalarField> impK() const;

            //- Return the implicit stiffness surface field
            //  This is the diffusivity for the Laplacian term
            tmp<surfaceScalarField> impKf() const;

            //- Return material residual i.e. a measured of how convergence of
            //  the material model
            scalar residual();

        // Edit

            //- Calculate the stress volField
            void correct(volSymmTensorField& sigma);

            //- Calculate the stress surfaceField
            void correct(surfaceSymmTensorField& sigma);

            //- Update the gradient of displacement vol field
            void grad
            (
                const volVectorField& D,
                volTensorField& gradD
            );

            //- Update the gradient of displacement vol field using pointD
            void grad
            (
                const volVectorField& D,
                const pointVectorField& pointD,
                volTensorField& gradD
            );

            //- Update the gradient of displacement surface field using pointD
            void grad
            (
                const volVectorField& D,
                const pointVectorField& pointD,
                surfaceTensorField& gradDf
            );

            //- Interpolate D to pointD
            void interpolate
            (
                const volVectorField& D,
                pointVectorField& pointD,
                const bool useVolFieldSigma = true
            );

            //- Return explicit Rhie-Chow correction field to quell oscillations
            tmp<volVectorField> RhieChowCorrection
            (
                const volVectorField& D,
                const volTensorField& gradD
            ) const;

            //- Update total accumulated fields
            void updateTotalFields();

            //- Lookup a field from the base mesh and interpolate it the subMesh
            template<class Type>
            tmp< GeometricField<Type, fvPatchField, volMesh> >
            lookupBaseMeshVolField
            (
                const word& fieldName,
                const fvMesh& subMesh
            ) const;

            //- Lookup a field from the base mesh and interpolate it the subMesh
            template<class Type>
            tmp< GeometricField<Type, fvsPatchField, surfaceMesh> >
            lookupBaseMeshSurfaceField
            (
                const word& fieldName,
                const fvMesh& subMesh
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "mechanicalModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
