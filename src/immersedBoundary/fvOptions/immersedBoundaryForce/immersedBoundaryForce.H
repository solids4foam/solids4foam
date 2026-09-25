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
    large K blocks every covered cell. With the volumeFraction weighting, the
    rate in the partially covered cells is also limited to
    \f$ \lambda/(1 - \lambda) C (\nu/w^2 + |U_i|/w) \f$, where w is the cell
    width and C the surface rate coefficient (default 3), so that it does
    not depend on the time step. The results do not depend on K for K of 1e3
    or more, the number of pressure correctors or the time step.

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

        // Coefficient of the penalisation rate in partially covered cells
        // with the volumeFraction weighting, optional, default 3; 0 for a
        // rate of lambda/((1 - lambda)*deltaT)
        surfaceRateCoeff 3;

        // Kinematic viscosity for the surface rate, optional, default the
        // nu entry of transportProperties
        nu              1e-3;

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
#include "immersedMarkers.H"

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

        //- Coefficient C of the penalisation rate in the partially covered
        //  cells with the volumeFraction weighting: the rate is
        //  lambda/(1 - lambda)*C*(nu/w^2 + |Ui|/w), where w is the cell
        //  width, so that it does not depend on the time step. The rate is
        //  lambda/((1 - lambda)*deltaT) if C is not positive
        scalar surfaceRateCoeff_;

        //- Kinematic viscosity for the surface rate
        scalar nu_;

        //- Has the viscosity been set?
        bool nuSet_;

        //- Penalisation rate (1/s) in each cell
        scalarField kappa_;

        //- Coupling coefficient of the incremental method
        scalar couplingCoeff_;

        // Marker method

            //- Marker spacing, relative to the local cell width
            scalar markerSpacing_;

            //- Retraction of the forcing points into the bodies, relative
            //  to the local cell width
            scalar markerRetraction_;

            //- Number of forcing iterations after the momentum predictor
            label nMarkerIterations_;

            //- Marker velocity response used in the forcing iterations:
            //  rowSum or diagonal
            word markerResponse_;

            //- Markers of each body
            PtrList<immersedMarkers> markers_;

            //- Forcing (acceleration) at the markers of each body, and at
            //  the start of the time step
            List<vectorField> markerForces_;
            List<vectorField> markerForces0_;

            //- Time index and outer corrector of the last forcing iteration
            label markerTimeIndex_;
            label markerCorr_;

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

        //- Kinematic viscosity for the surface rate: the nu entry, else the
        //  nu entry of transportProperties
        scalar nu();

        //- Width of a cell: its smallest bounding box width in the solved
        //  directions
        scalar cellWidth(const label celli) const;

        //- Move the bodies to the current time and calculate the occupancy
        void updateBodies(const bool move);

        //- Calculate the force on each body
        void calcForces(const volVectorField& U);

        //- Write the forces of the last completed time step
        void writeForces();

        //- Set the penalisation rate and the velocity of the bodies
        void setPenalty(const volVectorField& U);

        //- Multi-direct forcing iterations of the marker method, updating
        //  U, the forcing and the stored momentum matrix
        void markerForcing(volVectorField& U, const label nIter);


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
