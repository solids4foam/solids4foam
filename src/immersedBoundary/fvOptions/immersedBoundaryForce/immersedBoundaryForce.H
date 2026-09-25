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
    Foam::fv::immersedBoundaryForce

Description
    Hybrid fictitious domain-immersed boundary (HFDIB) method of openHFDIB
    and openHFDIB-DEM for bodies with a prescribed motion immersed in an
    incompressible flow, as a finite volume option for the momentum equation.

    Each body is a closed surface (see immersedBody) that covers some cells of
    the fluid mesh with an occupancy lambda, and moves with the velocity Ui.
    Two forcing methods are available.

    penalty (the default): the implicit volume penalisation

    \f[
        f = \kappa (U_i - U)
    \f]

    is added to the right-hand side of the momentum equation, with the
    \f$ -\kappa U \f$ part in the matrix, so that it also enters the pressure
    equation through 1/A. The penalisation rate is \f$ \kappa = K/\Delta t \f$
    in the fully covered cells, where K is the penalty coefficient. With the
    volumeFraction weighting (the default), \f$ \kappa \Delta t =
    \min(K, \lambda/(1 - \lambda)) \f$ in the partially covered cells, so
    that the penalised velocity is \f$ \lambda U_i + (1 - \lambda) U \f$; with
    the occupancy weighting, \f$ \kappa = K \lambda/\Delta t \f$, which for
    large K blocks every covered cell. The results do not depend on K for
    K of 1e3 or more, nor on the number of pressure correctors. With the
    occupancy weighting they do not depend on the time step either; with the
    volumeFraction weighting, the partially covered cells are relaxed towards
    Ui once per time step, so the results depend on the time step, and are
    closer to the body-fitted results for smaller time steps.

    incremental: the direct forcing of the openHFDIB-DEM pimpleHFDIBFoam
    solver. The explicit forcing f is added to the right-hand side of the
    momentum equation, and after each pressure corrector

    \f[
        f \leftarrow f + c \frac{U_i - U}{\Delta t}
    \f]

    in the covered cells, where \f$ c \f$ is the coupling coefficient; the
    updated forcing is used by the next pressure corrector. At the start of
    each time step, f is multiplied by the new occupancy; if the mesh moves
    again within the time step (moveMeshOuterCorrectors), f is set to zero in
    the cells no longer covered. The results depend on c, the time step and
    the number of pressure correctors. This method is kept for comparison, and
    reproduces pimpleHFDIBFoam with the vertexFraction occupancy.

    The force on each body is written every time step to
    postProcessing/<name>/<startTime>/<body>.dat, as the force -rho*sum(f*V)
    and torque exerted by the forcing, which exclude the inertia of the fluid
    inside the body, and that inertia rho*d/dt(sum(lambda*U*V)), whose sum
    with the force is the hydrodynamic force on a moving body (Uhlmann, 2005).

    The option requires a PIMPLE or PISO solution control on the mesh (e.g.
    pimpleFoam, or the solids4foam pimpleFluid fluid model), and is only
    available for incompressible flow.

    References:
    \verbatim
        Municchi, F. openHFDIB. DOI: 10.5281/zenodo.3940236.

        Isoz, M., Kotouč Šourek, M., Studeník, O., Kočí, P. (2022). Hybrid
        fictitious domain-immersed boundary solver coupled with discrete
        element method for simulations of flows laden with arbitrarily-shaped
        particles. Computers & Fluids, 244, 105538.
        DOI: 10.1016/j.compfluid.2022.105538.

        Studeník, O., Isoz, M., Kotouč Šourek, M., Kočí, P. (2024).
        OpenHFDIB-DEM: An extension to OpenFOAM for CFD-DEM simulations with
        arbitrary particle shapes. SoftwareX, 27.
        DOI: 10.1016/j.softx.2024.101871.

        Uhlmann, M. (2005). An immersed boundary method with direct forcing
        for the simulation of particulate flows. Journal of Computational
        Physics, 209, 448-476.
    \endverbatim

    Usage, in constant/fvOptions (and "libs (immersedBoundary);" in
    system/controlDict):
    \verbatim
    immersedBoundary
    {
        type            immersedBoundaryForce;

        // Forcing method, optional, default penalty
        method          penalty;

        // Penalty coefficient K, optional, default 1e3
        penaltyCoeff    1e3;

        // Weighting of the penalty in partially covered cells, optional,
        // default volumeFraction
        weighting       volumeFraction;

        // Coupling coefficient c of the incremental method, optional,
        // default 0.8
        couplingCoeff   0.8;

        // Occupancy threshold for a covered cell, optional, default 1e-4
        surfaceThreshold 1e-4;

        // Fluid density for the forces, optional, default the rho entry of
        // transportProperties, else 1 (kinematic forces)
        rho             1000;

        // Write the surfaces of moving bodies at write times, optional,
        // default yes
        writeSurfaces   yes;

        bodies
        {
            cylinder
            {
                surface     "cylinder.stl";
            }
        }
    }
    \endverbatim

    The fields lambda, Ui and f are written with the name of the option as a
    prefix, e.g. immersedBoundary:lambda. With the incremental method, f is
    read on restart.

Author
    Derived from the pimpleHFDIBFoam solver of openHFDIB-DEM
    (https://github.com/techMathGroup/openHFDIB-DEM, commit b4ea321), itself
    derived from openHFDIB (https://github.com/fmuni/openHFDIB).
    openHFDIB-DEM contributors: Federico Municchi (2016), Martin Isoz
    (2019-*), Martin Kotouč Šourek (2019-*), Ondřej Studeník (2020-*).
    openHFDIB and openHFDIB-DEM are distributed under the GNU Lesser General
    Public License version 3, whose terms allow this derived code to be
    distributed under the GNU General Public License version 3.

    Rewritten as a finite volume option for prescribed motion in solids4foam
    by Philip Cardiff, UCD, from the version adapted by Sairam Pamulaparthi
    Venkata, UCD, in cardiacFoam.

SourceFiles
    immersedBoundaryForce.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryForce_H
#define immersedBoundaryForce_H

#include "fvOption.H"
#include "volFields.H"
#include "PtrList.H"
#include "OFstream.H"
#include "IOdictionary.H"
#include "immersedBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pimpleControl;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class immersedBoundaryForce Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryForce
:
    public fv::option
{
    // Private Data

        //- Immersed bodies
        PtrList<immersedBody> bodies_;

        //- Forcing method: penalty (implicit volume penalisation, the
        //  default) or incremental (the explicit direct forcing of
        //  openHFDIB-DEM, increased after each pressure corrector)
        word method_;

        //- Penalty coefficient K: the penalisation rate in a fully covered
        //  cell is K/deltaT
        scalar penaltyCoeff_;

        //- Weighting of the penalisation in the partially covered cells:
        //  volumeFraction (the default), where the rate times deltaT is
        //  min(K, lambda/(1 - lambda)), so that the penalised velocity is
        //  lambda*Ui + (1 - lambda)*U, or occupancy, where the rate is
        //  K*lambda/deltaT
        word weighting_;

        //- Reconstruction of the target velocity in the cells next to the
        //  surface: none (the body velocity) or linear (a linear profile
        //  along the surface normal, through the body velocity on the
        //  surface and the velocity of a donor fluid cell)
        word reconstruction_;

        //- Velocity at the image point: cellPoint (interpolated, so that the
        //  target varies continuously as the body moves) or cell (the donor
        //  cell centre value)
        word imageInterpolation_;

        //- Penalisation rate (1/s) in each cell
        scalarField kappa_;

        //- If positive, the coefficient C of the penalisation rate in the
        //  cells next to the surface, C*(nu/w^2 + |Ub|/w) with w the cell
        //  width, which does not depend on the time step
        scalar bandRateCoeff_;

        //- Distance of the first image point from the surface, relative to
        //  the cell width
        scalar imageDistance_;

        //- Physical penalisation rate of a cell next to the surface
        scalar bandRate(const label celli, const vector& Ub) const;

        // Ghost cell method: the cells whose centre is inside a body are
        // penalised; the target velocity of those next to the fluid (the
        // ghost cells) is extrapolated linearly from a donor fluid cell
        // along the surface normal, so that the fluid velocity equals the
        // body velocity on the surface

            //- Penalised (inside) cells
            labelList insideCells_;

            //- Cells penalised in the previous configuration, which hold the
            //  body velocity rather than a fluid velocity, and so are not
            //  donors
            boolList penalisedOld_;

            //- Ghost cells
            labelList ghostCells_;

            //- Donor fluid cell of each ghost cell, -1 if none
            labelList ghostDonors_;

            //- Nearest surface point, outward unit normal, body velocity at
            //  the surface point and distance to it, for each ghost cell
            pointField ghostWallPoints_;
            vectorField ghostNormals_;
            vectorField ghostWallVelocities_;
            scalarField ghostDistances_;

            //- Ghost cells whose donor is on another processor: the ghost
            //  index, and the donor processor, for each
            labelList remoteGhosts_;
            labelList remoteDonorProcs_;

            //- Image point of each ghost cell with a local donor
            pointField ghostImagePoints_;

            //- Local donor cells sent to each processor, and their image
            //  points
            List<labelList> sendDonors_;
            List<pointField> sendImagePoints_;

            //- Remote ghost indices received from each processor, in the
            //  order they are sent
            List<labelList> recvGhosts_;

            //- Velocity and centre of the donor cell of each ghost cell,
            //  gathered from the donor processor
            vectorField ghostDonorU_;
            pointField ghostDonorC_;

        //- Coupling coefficient of the incremental method
        scalar couplingCoeff_;

        //- Occupancy threshold for a covered cell
        scalar surfaceThreshold_;

        //- Fluid density for the forces
        scalar rho_;

        //- Has the density been set?
        bool rhoSet_;

        //- Write the surfaces of moving bodies at write times
        bool writeSurfaces_;

        //- Occupancy (solid volume fraction)
        volScalarField lambda_;

        //- Velocity of the bodies in the covered cells, U elsewhere
        volVectorField Ui_;

        //- Forcing (acceleration)
        volVectorField f_;

        //- Forcing at the start of the time step
        vectorField f0_;

        //- Time index of the current configuration of the bodies
        label timeIndex_;

        //- Records the mesh points for which lambda was calculated
        IOdictionary pointsStamp_;

        //- Momentum matrix of the current outer corrector, set in constrain
        //  so that the pressure correctors see the updated forcing
        fvMatrix<vector>* UEqnPtr_;

        //- Time index and outer corrector of UEqnPtr_
        label UEqnTimeIndex_;
        label UEqnCorr_;

        //- Fluid momentum per unit density inside each body when the forces
        //  were last written, and the time index
        List<vector> momentum_;
        label momentumTimeIndex_;

        //- Fluid momentum per unit density inside each body at the end of
        //  the previous time step, if available
        List<vector> momentum0_;
        bool momentum0Valid_;

        //- Force files, on the master
        PtrList<OFstream> forceFiles_;

        //- Force, torque and inertia of the fluid inside each body, for the
        //  latest solution of the current time step, to be written once the
        //  time step is complete
        List<vector> force_;
        List<vector> torque_;
        List<vector> inertia_;

        //- Force on each body from the surface stress, sampled at image
        //  points along the surface normal
        List<vector> surfaceForce_;

        //- Image point distance for the surface force, if positive, else
        //  1.5 times the mean cell size of the reconstructed cells
        scalar surfaceForceDistance_;
        word forceTimeName_;
        bool forcesPending_;

        //- Has the missing momentum matrix been reported?
        bool warnedMatrix_;


    // Private Member Functions

        //- The solution control of the momentum equation
        const pimpleControl& pimple() const;

        //- Fluid density for the forces: the rho entry, else the rho entry
        //  of transportProperties, else 1. The fluid model may construct
        //  transportProperties after the option, so it is looked up on use
        scalar rho();

        //- Move the bodies to the current time and calculate the occupancy
        void updateBodies(const bool move);

        //- Calculate the force on each body
        void calcForces(const volVectorField& U);

        //- Write the forces of the last completed time step
        void writeForces();

        //- Force on a body from the pressure and viscous stress on its
        //  surface, sampled at image points at the distance l along the
        //  surface normal, with the wall pressure extrapolated with the body
        //  acceleration
        vector surfaceForce
        (
            const immersedBody& body,
            const volVectorField& U,
            const scalar l
        ) const;

        //- Set the penalisation rate and the velocity of the bodies
        void setPenalty(const volVectorField& U);

        //- Find the ghost cells and their donors, after the bodies are
        //  updated
        void findGhostCells();

        //- Find a donor fluid cell for the image points of a ghost cell,
        //  given its surface point, normal and width; -1 if none
        label findDonor
        (
            const point& wallPoint,
            const vector& normal,
            const scalar width,
            const volScalarField& inside,
            point& imagePoint
        ) const;

        //- Gather the velocity and centre of the donor cells
        void gatherDonors(const volVectorField& U);

        //- Reconstruct the target velocity in the cells next to the surface
        void reconstructTargets(const volVectorField& U);


public:

    //- Runtime type information
    TypeName("immersedBoundaryForce");


    // Constructors

        //- Construct from explicit source name and mesh
        immersedBoundaryForce
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        immersedBoundaryForce(const immersedBoundaryForce&) = delete;

        //- No copy assignment
        void operator=(const immersedBoundaryForce&) = delete;


    //- Destructor, which writes the forces of the last time step
    virtual ~immersedBoundaryForce();


    // Member Functions

        // Access

            //- Immersed bodies
            const PtrList<immersedBody>& bodies() const
            {
                return bodies_;
            }

            //- Occupancy
            const volScalarField& lambda() const
            {
                return lambda_;
            }

            //- Forcing
            const volVectorField& f() const
            {
                return f_;
            }


        // Evaluation

            //- Add the forcing to the momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Not implemented for compressible flow
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Store the momentum matrix
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Update the forcing after a pressure corrector
            virtual void correct(volVectorField& U);


        // IO

            //- Read the coefficients
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
