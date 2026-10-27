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
    meshMotionTotalDispSolid

Description
    Quasi-static solid model for moving a mesh, e.g. the fluid mesh with the
    solids4foamSolidModel mesh motion solver, which it shares with the mesh's
    owner.

    The unknown is the total displacement D from the reference points, but
    the equations are discretised on the current mesh, x = X + D, as
    displacementLaplacian does, so no copy of the reference mesh is needed.
    With finiteStrain yes, the kinematics are those of the current
    configuration:

        Finv = I - gradD^T,  F = inv(Finv),  J = det(F)

    where gradD is the gradient on the current mesh, and the momentum
    equation, div(sigma) = 0 for the Cauchy stress sigma, uses the current
    face areas. Once the mesh is at x = X + D, this is the spatial form of the
    total Lagrangian equations, so the mesh for a given boundary motion does
    not depend on how the boundary got there, for finite-strain laws too.
    With finiteStrain no (the default), the small-strain law is evaluated
    from gradD on the current mesh.

    The model has no inertia and no gravity, and it never updates the mesh:
    the mesh is moved by its owner, to the reference points plus pointD.
    Only history-independent (elastic) laws are supported.

    The equations are solved with the segregated corrector loop (the
    default, implicitSegregated) or, with solutionAlgorithm PETScSNES, with
    the PETSc nonlinear solver (SNES), whose residual is the same
    div(sigma) on the current mesh and whose preconditioner is the implicit
    stiffness Laplacian of the segregated loop. The PETSc options are read
    from the D entry of the fvSolution solvers of the moving mesh. The
    quasi-monolithic fluid-solid coupling requires PETScSNES.

    \verbatim
    solidModel      meshMotionTotalDisplacement;

    meshMotionTotalDisplacementCoeffs
    {
        // Finite-strain kinematics, for a finite-strain law, e.g.
        // diffusionHyperElastic (default no)
        finiteStrain        yes;

        nCorrectors         1000;
        rTol                1e-4;
        infoFrequency       100;

        // Optional: solve with PETSc SNES (default implicitSegregated)
        // solutionAlgorithm   PETScSNES;
        // optionsFile         petscOptions;
        // stopOnPetscError    yes;
    }
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    meshMotionTotalDispSolid.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_COM

#ifndef meshMotionTotalDispSolid_H
#define meshMotionTotalDispSolid_H

#include "solidModel.H"
#include "mechanicalConstitutiveLawManager.H"
#include "foamPetscSnesHelper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                   Class meshMotionTotalDispSolid Declaration
\*---------------------------------------------------------------------------*/

class meshMotionTotalDispSolid
:
    public solidModel,
    public foamPetscSnesHelper
{
    // Private data

        //- Use finite-strain kinematics
        const Switch finiteStrain_;

        //- Total deformation gradient
        volTensorField F_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Jacobian of the total deformation gradient
        volScalarField J_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the
        //  answer. It is updated at each solution, as the stiffness of the
        //  law can change with the mesh, e.g. through the stiffness scale
        //  factor of the diffusion laws
        volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;

        //- Reciprocal of impK
        volScalarField rImpK_;

        //- Number of unknowns per cell, used by the PETSc SNES solver
        const label blockSize_;


    // Private Member Functions

        //- Update gradD, and F, Finv and J, on the current mesh
        void updateKinematics();

        //- Update the cell-centred stress
        void correctStress();

        //- Update impK, impKf and rImpK from the law
        void updateImpK();

        //- Return the force at the faces of the current mesh, from sigma
        //  plus the stabilisation, with the traction conditions enforced;
        //  its divergence is the residual of the momentum equation
        tmp<surfaceVectorField> faceForce();

        //- Update the fields that follow from the solved D
        void updateDerivedFields();

        //- Evolve with the segregated corrector loop
        bool evolveImplicitSegregated();

        //- Evolve with the PETSc nonlinear solver (SNES)
        //  The residual is that of the segregated loop, and the implicit
        //  stiffness Laplacian is the preconditioner; a matrix-free
        //  (snes_mf_operator) Newton-Krylov method works well with it
        bool evolveSnes();

        //- Enforce traction boundary conditions on the provided traction
        //  field
        void enforceTractionBoundaries
        (
            surfaceVectorField& traction,
            const surfaceVectorField& n
        ) const;

        //- Disallow default bitwise copy construct
        meshMotionTotalDispSolid(const meshMotionTotalDispSolid&);

        //- Disallow default bitwise assignment
        void operator=(const meshMotionTotalDispSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        //  The Cauchy stress acts on the current faces, as for linear
        //  geometry
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("meshMotionTotalDisplacement");

    // Constructors

        //- Construct from components
        meshMotionTotalDispSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~meshMotionTotalDispSolid()
        {}


    // Member Functions

        // Access

            //- This model solves for D
            virtual volVectorField& solutionD()
            {
                return D();
            }

        // Edit

            //- Solve the mathematical model on the current mesh
            virtual bool evolve();

#ifdef USE_PETSC
            //- Initialise the Jacobian matrix
            virtual label initialiseJacobian(Mat& jac);

            //- Initialise the solution vector
            virtual label initialiseSolution(Vec& x);

            //- Calculate the residual given the solution
            virtual label formResidual
            (
                Vec f,         // Residual
                const Vec x    // Solution
            );

            //- Form the approximate Jacobian, i.e. the preconditioner
            virtual label formJacobian
            (
                Mat jac,       // Jacobian
                const Vec x    // Solution
            );
#endif // USE_PETSC

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // OPENFOAM_COM

// ************************************************************************* //
