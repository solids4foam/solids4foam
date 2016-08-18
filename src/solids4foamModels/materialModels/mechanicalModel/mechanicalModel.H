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

    To consider: cellZones would be more natural than a material indicator field
    but we would have to make sure no cell is in more than one cellZone.

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
#include "fvMeshSubset.H"

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

        //- Name of the cellZones for each material
        mutable wordList cellZoneNames_;

        //- Sub-mesh for each material
        //  Not used if there is only one material
        mutable PtrList<fvMeshSubset> subMeshes_;

        //- Stress volField for each mechanical law
        mutable PtrList<volSymmTensorField> subMeshSigma_;

        //- Stress surfaceField for each mechanical law
        mutable PtrList<surfaceSymmTensorField> subMeshSigmaf_;

    // Private Member Functions

        //- Make material sub-meshes
        void makeSubMeshes() const;

        //- Return a const reference to the subMeshes
        const PtrList<fvMeshSubset>& subMeshes() const;

        //- Check the cell zones are defined for each material
        void checkCellZones() const;

        //- Calculate the subMesh sigma fields
        void calcSubMeshSigma() const;

        //- Return the subMesh sigma volFields
        PtrList<volSymmTensorField>& subMeshSigma();

        //- Calculate the subMesh sigmaf fields
        void calcSubMeshSigmaf() const;

        //- Return the subMesh sigma surfaceFields
        PtrList<surfaceSymmTensorField>& subMeshSigmaf();

        //- Map a volField from the subMesh to the base mesh
        template<class Type>
        void mapSubMeshVolField
        (
            const label meshI,
            const GeometricField<Type, fvPatchField, volMesh>& subMeshField,
            GeometricField<Type, fvPatchField, volMesh>& baseMeshField
        ) const;

        //- Map a surfaceField from the subMesh to the base mesh
        template<class Type>
        void mapSubMeshSurfaceField
        (
            const label meshI,
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
                subMeshField,
            GeometricField<Type, fvsPatchField, surfaceMesh>& baseMeshField
        ) const;

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

            //- Update the gradient of displacement
            void grad
            (
                volTensorField& gradD, const volVectorField& D
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
