/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    beamModel

Description
    Virtual base class for beam models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Seevani Bali, UCD.

SourceFiles
    beamModel.C
    newBeamModel.C

\*---------------------------------------------------------------------------*/

#ifndef beamModel_H
#define beamModel_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dynamicFvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "labelScalarVector.H"
#include "labelScalarTimeVaryingVector.H"
#include "interpolationTable.H"
#include "crossSectionModel.H"
#include "labelPair.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class beamContactModel;

/*---------------------------------------------------------------------------*\
                           Class beamModel Declaration
\*---------------------------------------------------------------------------*/

class beamModel
:
    public IOdictionary
{
    // Private data

        //- Reference to the master time object
        Time& runTime_;

        //- Mesh, set only when this model constructs its own mesh
        autoPtr<dynamicFvMesh> meshOwnerPtr_;

        //- Mesh used by this model: owned, or borrowed from the registry
        const dynamicFvMesh* meshPtr_;

        //- Local to global cell addressing
        labelList localToGlobalCellAddressing_;

        //- Global to local cell addressing per processor
        List<labelList> globalToLocalCellAddressing_;

        //- Local to global point addressing per processor and beam
        PtrList<labelListList> localToGlobalBeamPointsAddressing_;

        //- Global to local point addressing per beam
        List<labelList> globalToLocalBeamPointsAddressing_;

        //- Cell face shared with upper neighbour cell
        mutable labelList* upperNeiCellFacesPtr_;

        //- Cell face shared with lower neighbour cell
        mutable labelList* lowerNeiCellFacesPtr_;

        //- Flow properties dictionary
        dictionary beamProperties_;

        //- Start patch index
        labelList startPatchIndex_;

        //- End patch index
        labelList endPatchIndex_;

        //- Index of first cell for every beam
        labelList startCells_;

        //- Number of cell for every beam
        labelList nBeamCells_;

        //- Beam cross-section model
        mutable PtrList<crossSectionModel> crossSections_;

        //- Beam cross-section radius
        scalarField R_;

        //- Beam convection velocity
        // scalarField U_;

        //- Young modulus
        // dimensionedScalar E_;
        scalarList E_;

        //- Shear modulus
        // dimensionedScalar G_;
        scalarList G_;

        //- Beam Density
        // dimensionedScalar rho_;
        scalarList rho_;

        //- Fluid density on Morrison's Equation
        dimensionedScalar rhof_;

        //- Gravitational Accl. ptr
        autoPtr<uniformDimensionedVectorField> gPtr_;

        // //- Beam cross-section area
        // dimensionedScalar A_;

        // //- Beam cross-section constant
        // dimensionedScalar I_;

        // //- Beam cross-section constant
        // dimensionedScalar J_;

        // //- Beam cross-section constant
        // dimensionedScalar EI_;

        // //- Beam cross-section constant
        // dimensionedScalar GJ_;

        // //- Beam cross-section constant
        // dimensionedScalar EA_;

        // //- Beam cross-section constant
        // dimensionedScalar GA_;

        //- Segment length
        mutable volScalarField L_;

        //- Distributed force
        volVectorField q_;

        //- Distributed moment
        volVectorField m_;

        //- Point forces
        labelScalarTimeVaryingVectorList pointForces_;

        //- Outer iteration counter
        label iOuterCorr_;

        //- Artificial scaling factor for changing the inertia
        //- tensor I_rho in the angular momentum tensor
        scalar kCI_;

        //- Conical pullys
        // mutable PtrList<conicalPulley> conicalPulleys_;

        // //- Toroidal pullys
        // mutable PtrList<toroidalPulley> toroidalPulleys_;

        // //- Beam contact model
        // autoPtr<beamContactModel> contactPtr_;

        //- celtaT time series
        interpolationTable<scalar> deltaTseries_;


    // Private Member Functions

        //- Create upper and lower cell faces
        void makeUpperLowerNeiCellFaces() const;

        //- Disallow default bitwise copy construct
        beamModel(const beamModel&);

        //- Disallow default bitwise assignment
        void operator=(const beamModel&);

protected:

    // Protected member functions

        //- Multibeam indicator
        tmp<surfaceScalarField> indicator(const label bI) const;

        // //- Return flow properties dictionary
        // dictionary& beamProperties()
        // {
        //     return beamProperties_;
        // }

public:

    //- Runtime type information
    TypeName("beamModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            beamModel,
            dictionary,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );

    // Constructors

        //- Construct from components
        beamModel
        (
            const word& type,
            Time& runTime,
            const word& region // = dynamicFvMesh::defaultRegion
        );

    // Selectors

        //- Select constructed from time
        static autoPtr<beamModel> New
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~beamModel();


    // Member Functions

        // Return global indices of local beam (spline) points
        tmp<labelField> globalPointsIndices(const label bI) const;

        // Return local indices of global beam (spline) points
        tmp<labelField> localPointsIndices(const label bI) const;

        // Access

            //- Return mesh
            const dynamicFvMesh& mesh() const
            {
                return *meshPtr_;
            }

            //- Return number of beams
            label nBeams() const
            {
                return mesh().cellZones().size();
            }

            //- Return number of cells (segments) for beam
            label nCells(const label bI) const
            {
                return mesh().cellZones()[bI].size();
            }

            // //- Return list of conical pulleys
            // const PtrList<conicalPulley>& conicalPulleys() const
            // {
  	        // return conicalPulleys_;
            // }

            // //- Return list of conical pulleys
            // const PtrList<toroidalPulley>& toroidalPulleys() const
            // {
  	        // return toroidalPulleys_;
            // }

            //- Return time
            const Time& runTime() const
            {
                return runTime_;
            }

            //- Return flow properties dictionary
            const dictionary& beamProperties() const
            {
                return beamProperties_;
            }

            //- Return start patch index
            const label& startPatchIndex(label bI = 0) const
            {
                return startPatchIndex_[bI];
            }

            //- Return end patch index
            const label& endPatchIndex(label bI = 0) const
            {
                return endPatchIndex_[bI];
            }

            //- Return start cells of beams
            const labelList& startCells() const
            {
                return startCells_;
            }

            //- Local-to-global beam points addressing
            const PtrList<labelListList>&
            localToGlobalBeamPointsAddressing() const
            {
                return localToGlobalBeamPointsAddressing_;
            }

            //- Local-to-global beam points addressing
            const labelListList&
            globalToLocalBeamPointsAddressing() const
            {
                return globalToLocalBeamPointsAddressing_;
            }

            //- For cell index return to which beam cell belongs
            label whichBeam(const label cellIndex) const;

            //- For cell index return local beam segment index
            label whichSegment(const label cellIndex) const;

            //- For local beam segment return global cell index
            label whichCell(const label beamIndex, const label segIndex) const;

            //- For global cell index return local cell index (parallel runs)
            // label localCellIndex(const label globalCellIndex) const;

            //- For global cell index return local cell index (parallel runs)
            labelPair procLocalCellIndex(const label globalCellIndex) const;

            //- For local cell index return local cell index (parallel runs)
            label globalCellIndex(const label localCellIndex) const;

            //- Beam cross-section model
             const PtrList<crossSectionModel>& crossSections() const
             {
                 return crossSections_;
             }

            //- Return beam cross-section radius
            const scalar& R(label bI = 0) const
            {
                return R_[bI];
            }

            //- Return beam convection velocity
            // const scalar& U(label bI = 0) const
            // {
            //     return U_
            //     [bI];
            // }

            inline dimensionedScalar A(label bI = 0) const
            {
                if (crossSections_.size())
                {
                    return
                        dimensionedScalar
                        (
                            "A", dimArea, crossSections_[bI].A()
                        );
                }

                return dimensionedScalar("A", dimArea, M_PI*sqr(R(bI)));
            }

            inline dimensionedScalar I(label bI = 0) const
            {
                if (crossSections_.size())
                {
                   return
                       dimensionedScalar
                       (
                           "I", dimArea*dimArea, crossSections_[bI].Ixx()
                       );
                }

                return dimensionedScalar("I", dimArea*dimArea, M_PI*pow(R(bI), 4)/4);
            }

            inline dimensionedScalar Iyy(label bI = 0) const
            {
                if (crossSections_.size())
                {
                   return
                       dimensionedScalar
                       (
                           "Iyy", dimArea*dimArea, crossSections_[bI].Ixx()
                       );
                }

                return dimensionedScalar("Iyy", dimArea*dimArea, M_PI*pow(R(bI), 4)/4);
            }

            inline dimensionedScalar Izz(label bI = 0) const
            {
                if (crossSections_.size())
                {
                   return
                       dimensionedScalar
                       (
                           "Izz", dimArea*dimArea, crossSections_[bI].Iyy()
                       );
                }

                return dimensionedScalar("Izz", dimArea*dimArea, M_PI*pow(R(bI), 4)/4);
            }

            inline dimensionedScalar J(label bI = 0) const
            {
                if (crossSections_.size())
                {
                   return
                       dimensionedScalar
                       (
                           "J", dimArea*dimArea, crossSections_[bI].IT()
                       );
                }

                return dimensionedScalar("J", dimArea*dimArea, M_PI*pow(R(bI), 4)/2);
            }

            const dimensionedScalar rho(label bI = 0) const
            {
                return dimensionedScalar("rhoBeam" + Foam::name(bI), dimDensity, rho_[bI]);
            }

            const dimensionedScalar& rhoFluid() const
            {
                return rhof_;
            }

            //- Return const reference to the gravitational acceleration
            UniformDimensionedField<vector>& g()
            {
                if (!gPtr_.valid())
                {
                    FatalError
                        << "g pointer for gravity not set" << abort(FatalError);
                }

                return gPtr_();
            }

            // const dimensionedScalar& E() const
            // {
            //     return E_;
            // }
            const dimensionedScalar E(label bI = 0) const
            {
                return dimensionedScalar("EBeam" + Foam::name(bI), dimPressure, E_[bI]);
            }

            // const dimensionedScalar& G() const
            // {
            //     return G_;
            // }

            const dimensionedScalar G(label bI = 0) const
            {
                return dimensionedScalar("GBeam" + Foam::name(bI), dimPressure, G_[bI]);
            }

            inline dimensionedScalar GA(label bI = 0) const
            {
                return G(bI)*A(bI);
            }

            inline dimensionedScalar GJ(label bI = 0) const
            {
                return G(bI)*J(bI);
            }

            inline dimensionedScalar EA(label bI = 0) const
            {
                return E(bI)*A(bI);
            }

            inline dimensionedScalar EI(label bI = 0) const
            {
                return E(bI)*I(bI);
            }

            inline dimensionedScalar EIyy(label bI = 0) const
            {
                return E(bI)*Iyy(bI);
            }

            inline dimensionedScalar EIzz(label bI = 0) const
            {
                return E(bI)*Izz(bI);
            }
            // Additional member functions added fo mass and rotation
 	          // inertia tensors for dynamics calculation
            inline dimensionedScalar ARho(label bI = 0) const
            {
                return A(bI)*rho(bI);
            }

            inline dimensionedScalar IyyRho(label bI = 0) const
            {
                return rho(bI)*Iyy(bI);
            }

 	        inline dimensionedScalar IzzRho(label bI = 0) const
            {
                return rho(bI)*Izz(bI);
            }

 	        inline dimensionedScalar JRho(label bI = 0) const
            {
                return rho(bI)*J(bI);
            }

            //- Return segment length
            const volScalarField& L() const
            {
                return L_;
            }

            //- Return segment length
            volScalarField& L()
            {
                return L_;
            }

            //- Return distributed force
            const volVectorField& q() const
            {
                return q_;
            }

            //- Return distributed force
            volVectorField& q()
            {
                return q_;
            }

            //- Return distributed force
            const volVectorField& m() const
            {
                return m_;
            }

            //- Return distributed force
            volVectorField& m()
            {
                return m_;
            }

            //- Return distributed force
            const labelScalarTimeVaryingVectorList& pointForces() const
            {
                return pointForces_;
            }

            //- Return outer iteration counter
            label& iOuterCorr()
            {
                return iOuterCorr_;
            }

            //- Return outer iteration counter
            const label& iOuterCorr() const
            {
                return iOuterCorr_;
            }

            //- Each beamModel must indicate whether W or DW is the primary
            //  solution variable
            virtual volVectorField& solutionW() const = 0;

            //- Each beamModel must indicate whether Theta or DTheta
            //  is the primary solution variable
            virtual volVectorField& solutionTheta() const = 0;

            //- Each beamModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix variable
            virtual surfaceTensorField& solutionLambda() const = 0;

            //- Each beamModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix variable
            virtual const volTensorField& solutionRMC() const = 0;

            // //- Each beamModel must indicate whether Lambda or DLambda
            // virtual const volVectorField& solutionCellTangent() const = 0;

            //- Return current beam points (ordered from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentBeamPoints
            (
                const label bI = 0
            ) const = 0;

            //- Return current beam points (ordered from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentGlobalBeamPoints
            (
                const label bI = 0
            ) const = 0;

            //- Return current beam points (ordered from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentBeamTangents
            (
                const label bI = 0
            ) const
            {
                // Not implemented
                return tmp<vectorField>(new vectorField(0, vector::zero));
            }

            //- Return current beam points (ordered from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentGlobalBeamTangents
            (
                const label bI = 0
            ) const
            {
                // Not implemented
                return tmp<vectorField>(new vectorField(0, vector::zero));
            }

            //- Return beam points and tangents
            virtual void currentGlobalBeamPointsAndTangents
            (
                const label bI,
                vectorField& points,
                vectorField& tangents
            ) const = 0;

            //- Return beam points and tangents
            virtual void currentGlobalBeamPointsAndTangents
            (
                const label bI,
                const labelListList& procSendPoints,
                // const labelListList& procReceivePoints,
                vectorField& points,
                vectorField& tangents
            ) const = 0;

            //- Final beam points
            tmp<vectorField> points(const label bI = 0) const;

            //- Return contact model
            // beamContactModel& contact();

            // //- Return contact model
            // const beamContactModel& contact() const;

            //- Check if beam-to-beam contact is active
            // bool contactActive() const
            // {
	        // return
            //     (
            //         beamContactActive()
            //      || pulleyContactActive()
            //     );
            // }

            // //- Check if pulley-to-beam contact is active
            // bool pulleyContactActive() const
            // {
	        // return
            //     (
            //         conicalPulleys().size()
            //      || toroidalPulleys().size()
            //     );
            // }

            // //- Check if beam-to-beam contact is active
            // bool beamContactActive() const
            // {
	        // return (mesh().cellZones().size() > 1);
            // }

            //- Extract beam point data from surface geometric field
            template<class Type>
            tmp<Field<Type> > beamPointData
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& sf,
                const label bI = 0
            ) const;

            //- Extract beam point data from surface geometric field
            template<class Type>
            tmp<Field<Type> > beamPointData
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& sf,
                const label bI,
		const label segI
            ) const;

            //- Extract segment data from volume geometric field
            template<class Type>
            tmp<Field<Type> > beamSegmentData
            (
                const GeometricField<Type, fvPatchField, volMesh>& vf,
                const label bI,
	        const label segI
            ) const;

            //- Return current time-step size
            virtual scalar deltaT() const;

        // Edit

            //- Evolve the solid solver
            virtual scalar evolve() = 0;

            //- Update total fields
            virtual void updateTotalFields()
            {
                // if (contactActive())
                // {
                //     contact().finalUpdate();
                // }
            }

            //- Update fields
            virtual void updateFields()
            {}

            //- Read fluidProperties dictionary
            virtual bool read();

            //- Write beam geometry in VTK format
            void writeVTK() const;

            //- Write fields
            virtual void writeFields();

            //- Artificial scaling factor for changing the inertia
            //- tensor I_rho in the angular momentum tensor
            scalar kCI() const
            {
                return kCI_;
            }

            //- Get upper neighbour cell faces
            const labelList& upperNeiCellFaces() const;

            //- Get lower neighbour cell faces
            const labelList& lowerNeiCellFaces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "beamModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
