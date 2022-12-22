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
    solidModel

Description
    Virtual base class for solid mechanics models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    solidModel.C
    solidModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solidModel_H
#define solidModel_H

#include "physicsModel.H"
#include "dynamicFvMesh.H"
#include "IOdictionary.H"
#include "scalarIOField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "thermalModel.H"
#include "mechanicalModel.H"
#include "nonLinearGeometry.H"
#include "uniformDimensionedFields.H"
#include "momentumStabilisation.H"
#include "globalPolyPatch.H"
#include "setCellDisplacements.H"
#include "OFstream.H"
#include "fvMatricesFwd.H"
#include "PackedBoolList.H"
#include "dualMeshToMeshMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                           Class solidModel Declaration
\*---------------------------------------------------------------------------*/

class solidModel
:
    public physicsModel,
    public regIOobject
{
    // Private data

        //- Mesh
        autoPtr<dynamicFvMesh> meshPtr_;

        //- Dual mesh
        mutable autoPtr<fvMesh> dualMeshPtr_;

        //- Dual mesh to primary mesh map
        mutable autoPtr<dualMeshToMeshMap> dualMeshToMeshMapPtr_;

        //- Solid properties dictionary
        IOdictionary solidProperties_;

        //- Derived solidModel type
        const word type_;

        //- Thermal model
        mutable autoPtr<thermalModel> thermalPtr_;

        //- Mechanical model
        mutable autoPtr<mechanicalModel> mechanicalPtr_;

        //- Total displacement field header
        IOobject Dheader_;

        //- Increment of displacement field header
        IOobject DDheader_;

        //- Point total displacement field header
        IOobject pointDheader_;

        //- Total displacement field
        volVectorField D_;

        //- Increment of displacement field
        //  DD = D - D.oldTime()
        volVectorField DD_;

        //- Velocity field
        volVectorField U_;

        //- Point mesh
        const pointMesh& pMesh_;

        //- Point total displacement field
        pointVectorField pointD_;

        //- Point increment of displacement field
        //  pointDD = pointD - pointD.oldTime()
        pointVectorField pointDD_;

        //- Gradient of total displacement
        volTensorField gradD_;

        //- Gradient of the displacement increment
        volTensorField gradDD_;

        //- Stress field
        //  This is the engineering stress for linear geometry approaches, and
        //  true (Cauchy) stress for nonlinear geometry (large strain)
        //  approaches
        volSymmTensorField sigma_;

        //- Density
        mutable autoPtr<volScalarField> rhoPtr_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField g_;

        //- Stabilisation term for momentum equation
        autoPtr<momentumStabilisation> stabilisationPtr_;

        //- Solution standard tolerance
        const scalar solutionTol_;

        //- Solution tighter tolerance
        const scalar alternativeTol_;

        //- Material law tolerance
        const scalar materialTol_;

        //- Write frequency for residuals information
        const int infoFrequency_;

        //- Maximum number of momentum correctors
        const int nCorr_;

        //- Minium number of momentum correctors
        const int minCorr_;

        //- Number of times the maximum number of correctors was reached in the
        //  momentum equation
        int maxIterReached_;

        //- Residual file
        autoPtr<OFstream> residualFilePtr_;

        //- Switch to write the residual field
        const Switch writeResidualField_;

        //- Switch to enable/disable enforce linear to help convergence
        Switch enforceLinear_;

        //- Under-relaxation method
        //  It can be 'fixed', 'Aitken' or 'QuasiNewton'
        const word relaxationMethod_;

        //- Aitken's under-relaxation factor field
        volScalarField aitkenAlpha_;

        //- Aitken's residual field
        volVectorField aitkenResidual_;

        //- Restart frequency for the Quasi-Newton method
        const int QuasiNewtonRestartFreq_;

        //- Quasi-Newton method input vector field
        DynamicList<vectorField> QuasiNewtonV_;

        //- Quasi-Newton method output vector field
        DynamicList<vectorField> QuasiNewtonW_;

        //- Quasi-Newton method time indices for the input-ouput fields
        DynamicList<scalar> QuasiNewtonT_;

        //- Reference D field used for Quasi-Newton method
        volVectorField DRef_;

        //- Reference D (without relaxation) field used for Quasi-Newton method
        volVectorField unrelaxedDRef_;

        //- Ptr list of global poly patches
        //  This is a copy of a patch which is reconstructed in full on each
        //  processor. This is convenient for interpolating between two patches
        //  in parallel, as used, for example, in the fluid-solid interaction
        //  procedure
        mutable PtrList<globalPolyPatch> globalPatchesPtrList_;

        //- Optional: set displacement at internal cells
        mutable autoPtr<setCellDisplacements> setCellDispsPtr_;

        //- Optional: restart flag
        //  Write out all fields required for a consistent restart
        const Switch restart_;

        //- Store rho*d2d2t(D) for efficiency
        mutable autoPtr<volVectorField> rhoD2dt2DPtr_;


    // Private Member Functions

        //- Disable solution in the out-of-plane direction for axisymmetric
        //  cases
        void checkWedges() const;

        //- Make the thermal model
        void makeThermalModel() const;

        //- Make the mechanical model
        void makeMechanicalModel() const;

        //- Make the rho field
        void makeRho() const;

        //- For calculating the position of the patch/face zone in the deformed
        //  configuration, we need to move the mesh by the point displacement
        //  field. For moving mesh approach (updated Lagrangian) we need the
        //  pointDD field, whereas for linear geometry (small strain) and total
        //  Lagrangian approaches, we need the pointD field. This function will
        //  check which field to return
        const pointVectorField& pointDorPointDD() const;

        //- Make setCellDisplacements
        void makeSetCellDisps() const;

        //- Return a reference to setCellDisplacements
        const setCellDisplacements& setCellDisps() const;

        //- Disallow default bitwise copy construct
        solidModel(const solidModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidModel&);


protected:

    // Protected member functions

        //- Make the dual mesh from the mesh
        void makeDualMesh() const;

        //- Mark features for dual mesh
        //  Copied from polyDualMeshApp.C
        void simpleMarkFeatures
        (
            const polyMesh& mesh,
            const PackedBoolList& isBoundaryEdge,
            const scalar featureAngle,
            const bool concaveMultiCells,
            const bool doNotPreserveFaceZones,
            labelList& featureFaces,
            labelList& featureEdges,
            labelList& singleCellFeaturePoints,
            labelList& multiCellFeaturePoints
        ) const;

        //- Return non-const reference to solid properties dictionary
        IOdictionary& solidProperties()
        {
            return solidProperties_;
        }

        //- Return non-const reference to solidModelCoeffs dictionary
        dictionary& solidModelDict();

        //- Return non-const reference to the thermal model
        thermalModel& thermal();

        //- Return non-const reference to the mechanical model
        mechanicalModel& mechanical();

        //- Return non-const reference to the density
        volScalarField& rho();

        //- Optional: cells can be forced to a specified displacement in the
        //  linear system
        void setCellDisps(fvVectorMatrix& DEqn);

        //- Apply fixed or Aitken's adaptive under-relaxation to the field
        void relaxField(volVectorField& D, int iCorr);

        //- Check if the equation has converged
        template<class Type>
        bool converged
        (
            const int iCorr,
            const scalar solverPerfInitRes,
            const int solverPerfNIters,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const bool writeResiduals = true
        );

        //- Return solution standard tolerance
        scalar solutionTol() const
        {
            return solutionTol_;
        }

        //- Return solution tighter tolerance
        scalar alternativeTol() const
        {
            return alternativeTol_;
        }

        //- Return material law tolerance
        scalar materialTol() const
        {
            return materialTol_;
        }

        //- Return write frequency for residuals information
        int infoFrequency() const
        {
            return infoFrequency_;
        }

        //- Return maximum number of correctors
        int nCorr() const
        {
            return nCorr_;
        }

        //- Return maxIterReached
        int& maxIterReached()
        {
            return maxIterReached_;
        }

        //- Return const reference to stabilisation
        const momentumStabilisation& stabilisation() const
        {
            return stabilisationPtr_();
        }

        //- Make rhoD2dt2D
        void makeRhoD2dt2D() const;

        //- Return non-const access to rhoD2dt2D
        volVectorField& rhoD2dt2D() const;


public:

    //- Runtime type information
    TypeName("solidModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidModel,
            dictionary,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );


    // Constructors

        //- Construct from components
        solidModel
        (
            const word& type,
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Selectors

        //- Select constructed from time
        static autoPtr<solidModel> New(Time& runTime, const word& region);


    // Destructor

        virtual ~solidModel();


    // Member Functions

        // Access

            //- Return mesh
            const dynamicFvMesh& mesh() const
            {
                return meshPtr_();
            }

            //- Return non-const mesh
            dynamicFvMesh& mesh()
            {
                return meshPtr_();
            }

            //- Return dualMesh
            const fvMesh& dualMesh() const;

            //- Return non-const dualMesh
            fvMesh& dualMesh();

            //- Return const dualMeshToMeshMap
            const dualMeshToMeshMap& dualMeshMap() const;

            //- Return point mesh
            const pointMesh& pMesh() const
            {
                return pMesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return meshPtr_().time();
            }

            //- Return solid properties dictionary
            const dictionary& solidProperties() const
            {
                return solidProperties_;
            }

            //- Return const-reference to the density
            const volScalarField& rho() const;

            //- Return const reference to solidModelCoeffs dictionary
            const dictionary& solidModelDict() const;

            //- Return const reference to the thermal model
            const thermalModel& thermal() const;

            //- Return const reference to the mechanical model
            const mechanicalModel& mechanical() const;

            //- Return non-const reference to D
            virtual volVectorField& D()
            {
                return D_;
            }

            //- Return const reference to D
            virtual const volVectorField& D() const
            {
                return D_;
            }

            //- Return non-const reference to DD
            virtual volVectorField& DD()
            {
                return DD_;
            }

            //- Return const reference to DD
            virtual const volVectorField& DD() const
            {
                return DD_;
            }

            //- Return non-const reference to U
            virtual volVectorField& U()
            {
                return U_;
            }

            //- Return const reference to U
            virtual const volVectorField& U() const
            {
                return U_;
            }

            //- Return non-const reference to pointD
            virtual pointVectorField& pointD()
            {
                return pointD_;
            }

            //- Return const reference to pointD
            virtual const pointVectorField& pointD() const
            {
                return pointD_;
            }

            //- Return non-const reference to pointDD
            virtual pointVectorField& pointDD()
            {
                return pointDD_;
            }

            //- Return const reference to pointDD
            virtual const pointVectorField& pointDD() const
            {
                return pointDD_;
            }

            //- Return velocity at a point
            virtual vector pointU(const label pointID) const;

            //- Return non-const reference to gradD
            virtual volTensorField& gradD()
            {
                return gradD_;
            }

            //- Return const reference to gradD
            virtual const volTensorField& gradD() const
            {
                return gradD_;
            }

            //- Return non-const reference to gradDD
            virtual volTensorField& gradDD()
            {
                return gradDD_;
            }

            //- Return const reference to gradDD
            virtual const volTensorField& gradDD() const
            {
                return gradDD_;
            }

            //- Return non-const reference to sigma
            virtual volSymmTensorField& sigma()
            {
                return sigma_;
            }

            //- Return const reference to sigma
            virtual const volSymmTensorField& sigma() const
            {
                return sigma_;
            }

            //- Return nonlinear geometry enumerator
            virtual nonLinearGeometry::nonLinearType nonLinGeom() const = 0;

            //- Return const reference to the gravitational acceleration
            const uniformDimensionedVectorField& g() const
            {
                return g_;
            }

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD() = 0;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                const label interfaceI
            ) const;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementOld
            (
                const label interfaceI
            ) const;

            //- Face zone acceleration
            virtual tmp<vectorField> faceZoneAcceleration
            (
                const label interfaceI
            ) const;

            //- Return enforceLinear switch
            const Switch& enforceLinear() const
            {
                return enforceLinear_;
            }

            //- Return enforceLinear switch
            Switch& enforceLinear()
            {
                return enforceLinear_;
            }

            //- Return restart switch
            const Switch& restart() const
            {
                return restart_;
            }

            //- Does this model take an incremental approach?
            //  i.e. does it solve for DD instead of D?
            //  This defaults to false but can be overwritten
            virtual bool incremental() const
            {
                return false;
            }

            //- Does this model move the mesh?
            virtual bool movingMesh() const
            {
                return false;
            }

            //- This function will check that the D field was read from disk
            virtual void DisRequired();

            //- This function will check that the DD field was read from disk
            virtual void DDisRequired();

            //- This function will check that the pointD field was read from disk
            virtual void pointDisRequired();

            //- Make the global poly patches
            void makeGlobalPatches
            (
                const wordList& patchNames,
                const bool currentConfiguration = true
            ) const;

            //- Return const-access to the global poly patches
            const PtrList<globalPolyPatch>& globalPatches() const;

            //- Clear global poly patches
            void clearGlobalPatches() const;

        // Edit

            //- Set traction at specified patch
            virtual void setTraction
            (
                fvPatchVectorField& tractionPatch,
                const vectorField& traction
            );

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label interfaceI,
                const label patchID,
                const vectorField& faceZoneTraction
            );

            //- Check if solid model is diverging using the Jacobian field
            virtual Switch& checkEnforceLinear(const volScalarField& J);

            //- Check if solid model is diverging using the Jacobian field
            virtual Switch& checkEnforceLinear(const surfaceScalarField& J);

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime)
            {}

            //- Evolve the solid model
            virtual bool evolve() = 0;

            //- Traction boundary surface normal gradient
            //  Given the user-specified boundary patch traction, this function
            //  should return the surface normal gradient
            //  This function is used by traction-type boundary conditions to
            //  set the boundary gradient
            virtual tmp<Foam::vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const = 0;

            //- Update total accumulated fields
            virtual void updateTotalFields();

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- Return the desired new time-step size
            virtual scalar newDeltaT();

            //- Move the mesh to the deformed configuration
            virtual void moveMesh
            (
                const pointField& oldPoints,
                const volVectorField& DD,
                pointVectorField& pointDD
            );

            //- End of the program
            virtual void end();

            //- Write function must be defined for regIOobjects
            virtual bool writeData(Ostream&) const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "solidModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
