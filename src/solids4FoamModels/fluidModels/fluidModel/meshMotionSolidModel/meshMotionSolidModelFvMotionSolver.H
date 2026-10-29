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

    The solid model, which must be meshMotionTotalDisplacement, is solved on
    the moving mesh itself, for the total displacement from the reference
    points, points0, and the mesh points are the reference points plus the
    solid model point displacement, pointD. As with displacementLaplacian, the
    equations are discretised on the current mesh, so no copy of the mesh is
    needed and, because the unknown is the total displacement and not an
    increment, the mesh for a given interface position does not drift with
    the path taken to reach it.

    The solid model reads its solidProperties and mechanicalProperties from
    constant/<region>, its D and pointD fields from <time>/<region>, and its
    schemes and solver settings from system/<region> of the moving mesh. Its
    fields whose names other models also use, e.g. U, p and rho, are private
    to it (see solidModel::sharedMesh()).

    The optional diffusivity entry, as for displacementLaplacian, sets the
    stiffnessScaleFactor field, which the diffusionElastic and
    diffusionHyperElastic mechanical laws use to scale their stiffness: the
    face diffusivity, limited to between minFactor and maxFactor times its
    average, divided by its average, and averaged to the cells. It is
    recalculated on the current mesh before each solution:

    \verbatim
    motionSolver    solids4foamSolidModel;

    // Optional
    diffusivity     quadratic inverseDistance (interface);
    maxFactor       100;
    minFactor       0.1;
    writeStiffScaleFactor no;
    \endverbatim

    If the solveEquations flag is false, the solid model pointD field is
    expected to be set by external functions.

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

        //- Mesh motion diffusivity, if specified
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Stiffness scale factor, if a diffusivity is specified
        autoPtr<volScalarField> stiffnessScaleFactorPtr_;

        //- Solid model
        mutable autoPtr<solidModel> modelPtr_;

        //- Flag to indicate that the model should solve its equations
        //  If this is false, it is assumed that the modelPtr_->pointD() field
        //  is updated by external functions
        mutable bool solveEquations_;


    // Private Member Functions

        //- Check that the solid model region is the region of the moving
        //  mesh, and return it
        word regionName() const;

        //- Create the stiffness scale factor, if a diffusivity is specified
        autoPtr<volScalarField> makeStiffnessScaleFactor();

        //- Update the stiffness scale factor on the current mesh
        void updateStiffnessScaleFactor();

        //- Create the solid model, which must be meshMotionTotalDisplacement
        void makeModel();

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
