/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    fluidSolidInterface

Description
    Virtual base class for fluid-solid interface coupling.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    fluidSolidInterface.C
    newFluidSolidInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSolidInterface_H
#define fluidSolidInterface_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fluidModel.H"
#include "solidModel.H"
#include "dynamicFvMesh.H"
#include "extendedGgiInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class fluidSolidInterface
:
    public IOdictionary
{
    // Private data

        //- FSI properties dictionary
        dictionary fsiProperties_;

        //- Fluid mesh
        dynamicFvMesh& fluidMesh_;

        //- Flow solver
        autoPtr<fluidModel> fluid_;

        //- Solid mesh
        fvMesh& solidMesh_;

        //- Solid solver
        autoPtr<solidModel> solid_;

        //- Solid interface patch index
        label solidPatchIndex_;

        //- Solid interface face zone index
        label solidZoneIndex_;

        //- Fluid interface patch index
        label fluidPatchIndex_;

        //- Fluid interface face zone index
        label fluidZoneIndex_;

        //- Solid face zone current points
        mutable vectorField* currentSolidZonePointsPtr_;

        //- Solid face zone primitive patch in current configuration
        mutable PrimitivePatch<face, List, const pointField&>*
            currentSolidZonePatchPtr_;

        //- Ggi zone-to-zone interpolation
        mutable extendedGgiZoneInterpolation* ggiInterpolatorPtr_;

        //- Outer correction loop stoping tolerance
        scalar outerCorrTolerance_;

        //- Maximal number of outer correctors
        label nOuterCorr_;

        //- Fsi coupling on/off
        Switch coupled_;

        //- Interface deformation limit
        scalar interfaceDeformationLimit_;

        //- Fluid zone point displacement
        vectorField fluidZonePointsDispl_;

        //- Fluid zone ref. point displacement
        vectorField fluidZonePointsDisplRef_;

        //- Fluid zone previous point displacement
        vectorField fluidZonePointsDisplPrev_;

        //- Solid zone point displacement
        vectorField solidZonePointsDispl_;

        //- Solid zone ref. point displacement
        vectorField solidZonePointsDisplRef_;

        //- Solid zone pressure
        scalarField solidZonePressure_;

        //- Current FSI residual
        vectorField residual_;

        //- Previous FSI residual
        vectorField residualPrev_;

        //- Maximal resudual norm
        scalar maxResidualNorm_;

        //- Outer corrector
        label outerCorr_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;

        //- Accumulated fluid side interface displacement
        mutable vectorIOField* accumulatedFluidInterfaceDisplacementPtr_;

        //- Min edge length for interface points at fluid side
        mutable scalarField* minEdgeLengthPtr_;

    // Private Member Functions

        //- Calculate current solid zone points
        void calcCurrentSolidZonePoints() const;

        //- Calculate current solid zone primitive patch
        void calcCurrentSolidZonePatch() const;

        //- Calculate fluid to solid ggi interpolator
        void calcGgiInterpolator() const;

        //- Accumulated fluid interface displacement
        void calcAccumulatedFluidInterfaceDisplacement() const;

        //- Calculate minimal edge lengths
        void calcMinEdgeLength() const;

        //- Return accumulated interface displacement
        vectorIOField& accumulatedFluidInterfaceDisplacement();

        //- Return minimal edge length
        const scalarField& minEdgeLength() const;

        //- Disallow default bitwise copy construct
        fluidSolidInterface(const fluidSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fluidSolidInterface&);


protected:

    // Protected member functions



public:

    //- Runtime type information
    TypeName("fluidSolidInterface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidSolidInterface,
            dictionary,
            (
                dynamicFvMesh& fluidMesh, fvMesh& solidMesh
            ),
            (fluidMesh, solidMesh)
        );


    // Constructors

        //- Construct from components
        fluidSolidInterface
        (
            const word& type,
            dynamicFvMesh& fluidMesh,
            fvMesh& solidMesh
        );


    // Selectors

        //- Select constructed from fluid and solid meshes
        static autoPtr<fluidSolidInterface> New
        (
            dynamicFvMesh& fluidMesh,
            fvMesh& solidMesh
        );


    // Destructor

        virtual ~fluidSolidInterface();


    // Member Functions

        // Access

            //- Return fluid mesh
            const dynamicFvMesh& fluidMesh() const
            {
                return fluidMesh_;
            }

            //- Return solid mesh
            const fvMesh& solidMesh() const
            {
                return solidMesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return fluidMesh_.time();
            }

            //- Return sub-dictionary
            const dictionary& fsiProperties() const
            {
                return fsiProperties_;
            }

            //- Return fluid model
            const fluidModel& fluid() const
            {
                return fluid_();
            }

            //- Return fluid model
            fluidModel& fluid()
            {
                return fluid_();
            }

            //- Return solid model
            const solidModel& solid() const
            {
                return solid_();
            }

            //- Return solid model
            solidModel& solid()
            {
                return solid_();
            }

            //- Return solid patch index
            label solidPatchIndex() const
            {
                return solidPatchIndex_;
            }

            //- Return solid face zone index
            label solidZoneIndex() const
            {
                return solidZoneIndex_;
            }

            //- Return fluid patch index
            label fluidPatchIndex() const
            {
                return fluidPatchIndex_;
            }

            //- Return fluid face zone index
            label fluidZoneIndex() const
            {
                return fluidZoneIndex_;
            }

            //- Return current solid zone points
            const vectorField& currentSolidZonePoints() const;

            //- Return current solid zone patch
            const PrimitivePatch<face, List, const pointField&>&
            currentSolidZonePatch() const;

            //- Return fluid to solid interpolator
            const extendedGgiZoneInterpolation& ggiInterpolator() const;

            //- Return outer corrector loop tolerance
            scalar outerCorrTolerance() const
            {
                return outerCorrTolerance_;
            }

            //- Return max numter of outer correctors
            label nOuterCorr() const
            {
                return nOuterCorr_;
            }

            //- Is it fluid and solid coupled
            const Switch& coupled() const
            {
                return coupled_;
            }

            //- Is it fluid and solid coupled
            Switch& coupled()
            {
                return coupled_;
            }

            //- Return relaxation factor
            scalar interfaceDeformationLimit() const
            {
                return interfaceDeformationLimit_;
            }

            //- Return fluid zone point displacement
            vectorField& fluidZonePointsDispl()
            {
                return fluidZonePointsDispl_;
            }

            //- Return fluid zone ref. point displacement
            vectorField& fluidZonePointsDisplRef()
            {
                return fluidZonePointsDisplRef_;
            }

            //- Return fluid zone previous point displacement
            vectorField& fluidZonePointsDisplPrev()
            {
                return fluidZonePointsDisplPrev_;
            }

            //- Return solid zone point displacement
            vectorField& solidZonePointsDispl()
            {
                return solidZonePointsDispl_;
            }

            //- Return solid zone ref. point displacement
            vectorField& solidZonePointsDisplRef()
            {
                return solidZonePointsDisplRef_;
            }

            //- Return solid zone pressure
            scalarField& solidZonePressure()
            {
                return solidZonePressure_;
            }

            //- Return current fsi residual
            vectorField& residual()
            {
                return residual_;
            }

            //- Return previous fsi residual
            vectorField& residualPrev()
            {
                return residualPrev_;
            }

            //- Return current outer iteration
            label& outerCorr()
            {
                return outerCorr_;
            }

            //- Return current outer iteration
            const label& outerCorr() const
            {
                return outerCorr_;
            }

        // Edit

            //- Initialize fields
            virtual void initializeFields();

            //- Update interpolator
            void updateInterpolator();

            //- Calculate interface displacement
            virtual void updateDisplacement() = 0;

            //- Move fluid mesh
            void moveFluidMesh();

            //- Update interface force
            virtual void updateForce();

            //- Update interface force
            virtual scalar updateResidual();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
