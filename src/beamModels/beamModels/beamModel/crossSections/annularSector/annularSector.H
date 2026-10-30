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
    annularSector

Description
    annularSector cross-section model

SourceFiles
    annularSector.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef annularSector_H
#define annularSector_H

#include "crossSectionModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "faceList.H"
#include "quadrilateral.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class annularSector Declaration
\*---------------------------------------------------------------------------*/

class annularSector
:
    public crossSectionModel
{
    // Private data

        //- Cross-section model dictionary
        dictionary crossSectionModelDict_;

        //- inner radius of the sector
        scalar innerRadius_;

        //- outerRadius of the sector
        scalar outerRadius_;

        //- angle of the sector
        scalar sectorAngle_;

        //- Number of segments of the sector in radial direction
    	label nRadialSegmentsSector_;

    	//- Cross-section area
        scalar A_;

        //- Cross-section second moment of area
        scalar Ixx_;

        //- Cross-section second moment of area
        scalar Iyy_;

        //- Number of circumferential points
        label nCircumferentialPoints_;

        //- Local mesh points
        vectorField points_;

        //- Local mesh faces
        // faceList quadFaces_;
        List<quadrilateral> quadFaces_;
        // faceList faces_;

    // Private Member Functions

        //- Calculate cross-section area
        scalar area() const;

        //- Calculate second moment of area
        void calcSecondMomentOfArea(scalar& Ixx, scalar& Iyy);

        //- Disallow copy construct
        annularSector(const annularSector&);

        //- Disallow default bitwise assignment
        void operator=(const annularSector&);

public:

    //- Runtime type information
    TypeName("annularSector");


    // Constructors

        //- Construct from dictionary
        annularSector
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~annularSector()
        {}


    // Member Functions


        //- Returns radius for that theta value
        virtual scalar radius(const scalar theta) const
        {
            return innerRadius_;   //Verify this!
        }

        //- Return cross-section radius
        virtual scalar R() const
        {			     // What does this do!
            return (outerRadius_ - innerRadius_);     //verify this with ZT!!
        }


        //- Return cross-section area
        virtual scalar A() const
        {
            return A_;
        }

        //- Return cross-section moment of inertia
        virtual scalar Ixx() const
        {
            return Ixx_;
        }

        //- Return cross-section moment of inertia
        virtual scalar Iyy() const
        {
            return Iyy_;
        }

        //- Return cross-section Saint-Venant torsion constant
        virtual scalar IT() const
        {
            return Ixx() + Iyy();
        }

        // //- Return cross-section coefficient of proportionality
        // //  between fully plastic axial force and fully plastic moment
        // virtual scalar CNM() const
        // {
        //     return 4.0/innerRadius_;  //Not sure what this is!
        // }

        //- Circumferential points
        virtual tmp<vectorField> circumferentialPoints() const;

        //- Number of circumferential points
        virtual label nCircumferentialPoints() const
        {
            return nCircumferentialPoints_;
        }

        //- Structured mesh points
        virtual const vectorField& points() const
        {
            return points_;
        }

        // //- Return Green-Lagrangian strain at structured mesh points
        // //  for provided average tranlational and rotational strains
        // virtual tmp<vectorField> greenLagrangianStrain
        // (
        //     const vector& Gamma,
        //     const vector& K
        // ) const;

        // //- Return cross-section resultant force
        // virtual tmp<scalarField> resultantTotalForce
        // (
        //     const vectorField& S
        // ) const;

        // //- Return cross-section resultant force
        // virtual vector resultantForce(const vectorField& S) const;

        // //- Return cross-section resultant moment
        // virtual vector resultantMoment(const vectorField& S) const;

        // //- Return cross-section resultant tangent matrices
        // virtual void resultantTangentMatrices
        // (
        //     const diagTensorField& dSdE,
        //     tensor& DQDGamma,
        //     tensor& DQDK,
        //     tensor& DMDGamma,
        //     tensor& DMDK
        // ) const;

        // //- Write cross-section scalar field in VTK format
        // virtual void writeVTK
        // (
        //     const fileName& fn,
        //     const scalarField& sf0,
        //     const scalarField& sf1,
        //     const vectorField& vf
        // ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
