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
    rectangle

Description
    Rectangle cross-section model

SourceFiles
    rectangle.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef rectangle_H
#define rectangle_H

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
                         Class rectangle Declaration
\*---------------------------------------------------------------------------*/

class rectangle
:
    public crossSectionModel
{
    // Private data

        //- Cross-section model dictionary
        dictionary crossSectionModelDict_;

        //- Rectangle length
        scalar b_;

        //- Rectangle width
        scalar h_;

        //- Number of segments in x-direction
        scalar nx_;

        //- Number of segments in y-direction
        scalar ny_;

        //- Artificial area scaling factor which affects the mechanical
        //- properties of the beams
        scalar scalingArea_;

        //- Artificial second moment of area scaling factor which affects
        //- the mechanical properties of the beams
        scalar scalingSecondMomentArea_;

        //- Local mesh points
        vectorField points_;

        //- Local mesh faces
        // faceList quadFaces_;
        List<quadrilateral> quadFaces_;
        // faceList faces_;

    // Private Member Functions

        //- Disallow copy construct
        rectangle(const rectangle&);

        //- Disallow default bitwise assignment
        void operator=(const rectangle&);

public:

    //- Runtime type information
    TypeName("rectangle");


    // Constructors

        //- Construct from dictionary
        rectangle
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~rectangle()
        {}


    // Member Functions

        //- Return cross-section area
        virtual scalar radius(const scalar theta) const
        {
            return h_/2;
        }

        //- Return cross-section area
        virtual scalar R() const
        {
            return h_/2;
        }

        //- Return the width b
        virtual scalar b() const
        {
            return b_;
        }

        //- Return cross-section area
        virtual scalar A() const
        {
            return scalingArea_*b_*h_;
        }

        //- Return cross-section moment of inertia
        virtual scalar Ixx() const
        {
            return scalingSecondMomentArea_*b_*pow(h_, 3)/12;
        }

        //- Return cross-section moment of inertia
        virtual scalar Iyy() const
        {
            return scalingSecondMomentArea_*h_*pow(b_, 3)/12;
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
        //     return 4.0/h_;
        // }

        //- Circumferential points
        virtual tmp<vectorField> circumferentialPoints() const
        {
            label nTheta = 4;

            tmp<vectorField> tCircumPoints
            (
                new vectorField(nTheta, vector::zero)
            );

            // Counter-clock-wise orientation
            tCircumPoints.ref()[0] = vector(0, b_/2, -h_/2);
            tCircumPoints.ref()[1] = vector(0, -b_/2, -h_/2);
            tCircumPoints.ref()[2] = vector(0, -b_/2, h_/2);
            tCircumPoints.ref()[3] = vector(0, b_/2, h_/2);

            return tCircumPoints;
        }

        //- Number of circumferential points
        virtual label nCircumferentialPoints () const
        {
            label nTheta = 4;
            return nTheta;
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
        // void resultantTangentMatrices
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
