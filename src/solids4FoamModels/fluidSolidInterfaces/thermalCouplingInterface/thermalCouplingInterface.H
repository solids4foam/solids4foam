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
    thermalCouplingInterface

Description
    Strong Dirichlet-Neumann coupling with fixed under-relaxation.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    thermalCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef thermalCouplingInterface_H
#define thermalCouplingInterface_H

#include "fixedRelaxationCouplingInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                  Class thermalCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class thermalCouplingInterface
:
    public fixedRelaxationCouplingInterface
{
    // Private data

        //- Fixed under-relaxation factor
        const scalar relaxationFactor_;

        //- Predict temperature and heat-flux at fluid and solid interface
        const Switch predictTemperatureAndHeatFlux_;

        //- Flag to enable mechanical coupling between the fluid and solid
        const Switch mechanicalCoupling_;

        //- Residual file pointer
        autoPtr<OFstream> thermalResidualFilePtr_;

        //- List of maximal thermal residuals norm
        List<scalar> maxThermalResidualsNorm_;

        //- Old-time solid face-zone temperature
        mutable PtrList<scalarField> oldSolidFaceZoneTemperature_;

        //- Old-time solid face-zone heat flux
        mutable PtrList<scalarField> oldSolidFaceZoneHeatFlux_;

        //- Old-old-time solid face-zone temperature
        mutable PtrList<scalarField> oldOldSolidFaceZoneTemperature_;

        //- Old-old-time solid face-zone heat flux
        mutable PtrList<scalarField> oldOldSolidFaceZoneHeatFlux_;

        //- Time index
        label timeIndex_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalCouplingInterface
        (
            const thermalCouplingInterface&
        );

        //- Disallow default bitwise assignment
        void operator=(const thermalCouplingInterface&);


public:

    //- Runtime type information
    TypeName("thermal");


    // Constructors

        //- Construct from components
        thermalCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~thermalCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Update temperature and heat flux on fluid side
            void updateHeatFluxAndTemperatureOnFluidInterface();

            //- Update temperature and heat flux on solid side
            void updateHeatFluxAndTemperatureOnSolidInterface();

            //- Return reference to the residualFile
            OFstream& thermalResidualFile();

            //- Calculate thermal residual
            scalar calcThermalResidual();

            //- Return the current rate heat flow on the interface
            //  given the zone and zone traction field
            scalar heatFlowRateOnInterface
            (
                const standAlonePatch& zone, const scalarField& zoneTraction
            ) const;

            //- Set equivalent interface heat transfer coeff (delta/lambda_eq)
            void setEqInterHeatTransferCoeff();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
