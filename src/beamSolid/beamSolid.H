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

\*---------------------------------------------------------------------------*/

#ifndef beamSolid_H
#define beamSolid_H

#include "solidModel.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                         Class beamSolid Declaration
\*---------------------------------------------------------------------------*/

class beamSolid
:
    public solidModel
{
    // Private data

        //- Beam mechanics on the solid model's mesh
        autoPtr<beamModel> beam_;

    // Private Member Functions

        //- Synchronise the solid displacement interface with beam kinematics
        void updateDisplacements();

public:

    //- Runtime type information
    TypeName("beamSolid");

    // Constructors

        //- Construct from time and region
        beamSolid(Time& runTime, const word& region);

    // Destructor

        virtual ~beamSolid() = default;

    // Member Functions

        //- Total Lagrangian beam kinematics
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }

        //- Primary displacement unknown remains W
        virtual volVectorField& solutionD()
        {
            return beam_().solutionW();
        }

        //- Beam translational velocity
        virtual volVectorField& U();

        //- Beam translational velocity
        virtual const volVectorField& U() const;

        //- Solve beam translation and rotation
        virtual bool evolve();

        //- Update beam and solid-interface displacement fields
        virtual void updateTotalFields();

        //- Write beam fields without continuum stress post-processing
        virtual void writeFields(const Time& runTime);

        //- Continuum traction boundary conditions are not beam loads
        virtual tmp<vectorField> tractionBoundarySnGrad
        (
            const vectorField& traction,
            const scalarField& pressure,
            const fvPatch& patch
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

#endif

// ************************************************************************* //
