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

    The optional diffusivity entry, as for displacementLaplacian, creates the
    stiffnessScaleFactor field, which the diffusionElastic and
    diffusionHyperElastic mechanical laws use to scale their stiffness. It is
    calculated once, from the mesh at construction, as the face diffusivity
    limited to between minFactor and maxFactor times its average, divided by
    its average, and averaged to the cells:

    \verbatim
    motionSolver    solids4foamSolidModel;

    regionName      fluid;

    // Optional
    diffusivity     quadratic inverseDistance (interface);
    maxFactor       100;
    minFactor       0.1;
    writeStiffScaleFactor no;
    \endverbatim

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

        //- Stiffness scale factor for the diffusionElastic and
        //  diffusionHyperElastic mechanical laws, if a diffusivity is
        //  specified
        autoPtr<volScalarField> stiffnessScaleFactorPtr_;

        //- Solid model
        mutable autoPtr<solidModel> modelPtr_;

        //- Flag to indicate that the model should solve its equations
        //  If this is false, it is assumed that the modelPtr_->pointD() field
        //  is updated by external functions
        mutable bool solveEquations_;


    // Private Member Functions

        //- Calculate the stiffness scale factor, if a diffusivity is
        //  specified
        autoPtr<volScalarField> makeStiffnessScaleFactor() const;

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
