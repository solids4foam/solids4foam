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
    Foam::mechanicalModel

Description
    The mechanical model class takes care of reading the mechanical properties
    and creating mechanical laws.

    The solidSubMeshes class takes care of multiple mechanical laws, where each
    each material region has its own mechanical law.

    Corrections are applied at bi-material interfaces to ensure continuity of
    stress without oscillations.

SourceFiles
    mechanicalModel.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalModel_H
#define mechanicalModel_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "mechanicalLaw.H"
#include "Switch.H"
#include "solidSubMeshes.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "enhancedVolPointInterpolation.H"
#else
    #include "newLeastSquaresVolPointInterpolation.H"
#endif
#include "nonLinearGeometry.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "hashedWordList.H"

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

        //- Stress field
        volSymmTensorField sigma_;

        //- Is the current case plane stress or plain strain/3-D
        //  Todo: this should be moved to actual mechanical law
        const Switch planeStress_;

        //- Does the solidModel use an incremental approach?
        //  i.e. does the solidModel solve for DD instead of D
        const bool incremental_;

        //- Name of the cellZone for each material
        mutable wordList cellZoneNames_;

        //- Sub-meshes for each mechanicalLaw
        mutable autoPtr<solidSubMeshes> solSubMeshes_;

        //- List with the names of all internal state variables
        //  Key = datatypes (volScalarField, etc.)
        //  Values = hashed wordLists of names of internal vars
        mutable HashTable<hashedWordList> internalVariables_;

        //- Flag to see if internal variables got read
        bool gotInternalVariables_;

        //- Which datatypes will be considered internal state variables
        static const hashedWordList internalVariableTypes;

    // Private Member Functions

        //- Make material sub-meshes
        void makeSolSubMeshes() const;

        //- Get the internal variable names
        bool getInternalVariables();

        //- Clear out demand driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        mechanicalModel(const mechanicalModel&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalModel&);


public:

    // Constructors

        //- Construct from an fvMesh
        mechanicalModel
        (
            const fvMesh& mesh,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const bool incremental = false
        );


    //- Destructor
    ~mechanicalModel();


    // Member Functions

        // Access

            //- Return const access to the mesh
            const fvMesh& mesh() const;

            //- Return reference to the stress field
            volSymmTensorField& sigma();

            //- Return const reference to the stress field
            const volSymmTensorField& sigma() const;

            //- Return a const reference to the vol-to-point base mesh
            //  interpolator
#ifdef OPENFOAM_NOT_EXTEND
            const enhancedVolPointInterpolation& volToPoint() const;
#else
            const newLeastSquaresVolPointInterpolation& volToPoint() const;
#endif

            //- Return a const reference to the solidSubMeshes
            const solidSubMeshes& solSubMeshes() const;

            //- Return a const reference to the solidSubMeshes
            solidSubMeshes& solSubMeshes();

            //- Return initial density
            tmp<volScalarField> rho() const;

            //- Return the implicit stiffness
            //  This is the diffusivity for the Laplacian term
            tmp<volScalarField> impK() const;

            //- Return the implicit stiffness surface field
            //  This is the diffusivity for the Laplacian term
            tmp<surfaceScalarField> impKf() const;

            //- Return bulk modulus
            tmp<volScalarField> bulkModulus() const;

            //- Return material residual i.e. a measured of how convergence of
            //  the material model
            scalar residual();

        // Edit

            //- Calculate the stress volField
            void correct(volSymmTensorField& sigma);

            //- Calculate the stress surfaceField
            void correct(surfaceSymmTensorField& sigma);

            //- Calculate the stress pointField given the displacement gradient
            void correct
            (
                pointSymmTensorField& sigma, const pointTensorField& gradD
            );

            //- Map grad field to subMeshes
            void mapGradToSubMeshes(const volTensorField& gradD);

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

            //- Update the gradient of displacement vol and surface field using
            //  pointD
            void grad
            (
                const volVectorField& D,
                const pointVectorField& pointD,
                volTensorField& gradD,
                surfaceTensorField& gradDf
            );

            //- Interpolate D to pointD
            void interpolate
            (
                const volVectorField& D,
                pointVectorField& pointD,
                const bool useVolFieldSigma = true
            );

            //- Interpolate D to pointD and use gradD for
            //  enhancedVolPointInterpolation
            //  Note: for multi-material cases, it is assumed that gradD is
            //  already up-to-date on sub-meshes
            void interpolate
            (
                const volVectorField& D,
                const volTensorField& gradD,
                pointVectorField& pointD,
                const bool useVolFieldSigma = true
            );

            //- Update total accumulated fields
            void updateTotalFields();

            //- Return the desired new time-step size
            scalar newDeltaT();

            //- Move subMeshes
            void moveSubMeshes();

            //- Set restart: some mechanical laws may change the writeOpt of
            //  fields accordingly
            void setRestart();

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- Write current state of 'materialProperties' dictionary to disk
            void writeDict();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
