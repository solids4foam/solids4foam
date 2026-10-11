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
    Foam::meshMotionSolidModelFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh using a solids4foam solidModel.

    The solid model is solved on the region given by regionName. When this
    is a separate region whose mesh is a copy of the moving mesh, e.g.
    constant/fluidMeshMotion/polyMesh copied from constant/fluid/polyMesh,
    the solid model is solved on this static reference mesh for the total
    displacement from the reference points, and the mesh points are the
    reference points plus the solid model point displacement. The meshes must
    have the same points, in the same order, including in parallel, where
    they must be decomposed in the same way.

    When the reference mesh is used, the optional diffusivity entry, as for
    displacementLaplacian, creates the stiffnessScaleFactor field on the
    reference mesh, which the diffusionElastic and diffusionHyperElastic
    mechanical laws use to scale their stiffness; it is calculated once, from
    the reference mesh, as the face diffusivity limited to between minFactor
    and maxFactor times its average, divided by its average, and averaged to
    the cells:

    \verbatim
    motionSolver    solids4foamSolidModel;

    regionName      fluidMeshMotion;

    // Optional
    diffusivity     quadratic inverseDistance (interface);
    maxFactor       100;
    minFactor       0.1;
    writeStiffScaleFactor no;
    \endverbatim

    When regionName is the region of the moving mesh itself, the solid model
    shares the moving mesh; its pointD field is then expected to be set by
    external functions (see solveEquations).

SourceFiles
    meshMotionSolidModelFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef meshMotionSolidModelFvMotionSolver_H
#define meshMotionSolidModelFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class meshMotionSolidModelFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class meshMotionSolidModelFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Static reference mesh on which the solid model is solved, if the
        //  solid model region is not the region of the moving mesh
        autoPtr<dynamicFvMesh> referenceMeshPtr_;

        //- Stiffness scale factor on the reference mesh
        autoPtr<volScalarField> stiffnessScaleFactorPtr_;

        //- Solid model
        mutable autoPtr<solidModel> modelPtr_;

        //- Flag to indicate that the model should solve its equations
        //  If this is false, it is assumed that the modelPtr_->pointD() field
        //  is updated by external functions
        mutable bool solveEquations_;


    // Private Member Functions

        //- Read the reference mesh, unless the solid model region is the
        //  region of the moving mesh
        autoPtr<dynamicFvMesh> makeReferenceMesh() const;

        //- Calculate the stiffness scale factor on the reference mesh, if a
        //  diffusivity is specified
        autoPtr<volScalarField> makeStiffnessScaleFactor() const;

        //- Check that the reference mesh points are the reference points
        void checkReferenceMesh() const;

        //- No copy construct
        meshMotionSolidModelFvMotionSolver
        (
            const meshMotionSolidModelFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const meshMotionSolidModelFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("solids4foamSolidModel");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        meshMotionSolidModelFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        meshMotionSolidModelFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~meshMotionSolidModelFvMotionSolver();


    // Member Functions

        //- Access to the solid model
        const solidModel& model() const
        {
            return modelPtr_();
        }

        //- Non-const access to the solid model
        solidModel& model()
        {
            return modelPtr_();
        }

        //- Return true if the solid model is solved on a separate, static
        //  reference mesh
        bool referenceMesh() const
        {
            return referenceMeshPtr_.valid();
        }

        //- Access to the solveEquations flag
        bool& solveEquations() const
        {
            return solveEquations_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
