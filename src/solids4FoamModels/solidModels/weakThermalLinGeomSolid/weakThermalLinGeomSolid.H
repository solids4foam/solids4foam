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
    weakThermalLinGeomSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the heat equation is solved.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    weakThermalLinGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef weakThermalLinGeomSolid_H
#define weakThermalLinGeomSolid_H

#include "linGeomTotalDispSolid.H"
#include "thermalModel.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "SolverPerformance.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class weakThermalLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class weakThermalLinGeomSolid
:
    public linGeomTotalDispSolid
{
    // Private data

        //- Thermal model
        thermalModel thermal_;

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField k_;

        //- Temperature
        volScalarField T_;

        //- Temperature gradient
        volVectorField gradT_;

        //- Absolute temperature tolerance (in degrees)
        const scalar absTTol_;

    // Private Member Functions

        //- Check convergence
        bool converged
        (
            const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
            const SolverPerformance<scalar>& solverPerfT,
#else
            const lduSolverPerformance& solverPerfT,
#endif
            const volScalarField& T
        );

    //- Disallow default bitwise copy construct
        weakThermalLinGeomSolid(const weakThermalLinGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const weakThermalLinGeomSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("weakThermalLinearGeometry");

    // Constructors

        //- Construct from components
        weakThermalLinGeomSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~weakThermalLinGeomSolid()
        {}


    // Member Functions

        // Access

            //- Return const access to the cell temperature field
            const volScalarField& T() const
            {
                return T_;
            }

            //- Return non-const access to the cell temperature field
            volScalarField& T()
            {
                return T_;
            }

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
