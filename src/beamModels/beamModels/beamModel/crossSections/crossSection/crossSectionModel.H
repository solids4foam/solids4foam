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
    crossSectionModel

Description
    Cross-section model

SourceFiles
    crossSectionModel.C
    newCrossSectionModel.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef crossSectionModel_H
#define crossSectionModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "vectorField.H"
#include "tensorField.H"
#include "diagTensorField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------* \
                         Class crossSectionModel Declaration
\*---------------------------------------------------------------------------*/

class crossSectionModel
{
    // Private data

        //- Name
        const word name_;

    // Private Member Functions

        //- Disallow copy construct
        crossSectionModel(const crossSectionModel&);

        //- Disallow default bitwise assignment
        void operator=(const crossSectionModel&);

public:

    //- Runtime type information
    TypeName("crossSectionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            crossSectionModel,
            dictionary,
            (
                const word name,
                const dictionary& dict
            ),
            (
                name,
                dict
            )
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<crossSectionModel> New
        (
            const word& name,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        crossSectionModel
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~crossSectionModel()
        {}


    // Member Functions

        //- Tensor to scalar rectangular matrix
        void tensorToMatrix
        (
            scalarRectangularMatrix& m,
            const tensor& t
        ) const;

        //- Tensor to scalar rectangular matrix
        void tensorToMatrix
        (
            scalarRectangularMatrix& m,
            const diagTensor& t
        ) const;

        //- Calculate matrix A
        void calcA
        (
            scalarRectangularMatrix& A,
            const vector& p
        ) const;

        //- Return cross-section radius
        virtual scalar radius(const scalar theta) const = 0;

        //- Return polar angel for given point in cartesion coordinates
        scalar polarAngle(const scalar x, const scalar y) const;

        //- Return cross-section equivalent radius
        virtual scalar R() const = 0;

        //- Return cross-section area
        virtual scalar A() const = 0;

        //- Return cross-section moment of inertia
        virtual scalar Ixx() const = 0;

        //- Return cross-section moment of inertia
        virtual scalar Iyy() const = 0;

        //- Return cross-section Saint-Venant torsion constant
        virtual scalar IT() const = 0;

        // //- Return cross-section coefficient of proportionality
        // //  between fully plastic axial force and fully plastic moment
        // virtual scalar CNM() const = 0;

        //- Circumferential points
        virtual tmp<vectorField> circumferentialPoints() const = 0;

        //- Number of circumferential points
        virtual label nCircumferentialPoints () const = 0;

        //- Structured mesh points
        virtual const vectorField& points() const = 0;

        // //- Return Green-Lagrangian strain at structured mesh points
        // //  for provided average tranlational and rotational strains
        // virtual tmp<vectorField> greenLagrangianStrain
        // (
        //     const vector& Gamma,
        //     const vector& K
        // ) const = 0;

        // //- Return cross-section resultant force
        // virtual tmp<scalarField> resultantTotalForce
        // (const vectorField& S) const = 0;

        // //- Return cross-section resultant force
        // virtual vector resultantForce(const vectorField& S) const = 0;

        // //- Return cross-section resultant moment
        // virtual vector resultantMoment(const vectorField& S) const = 0;

        // //- Return cross-section resultant tangent matrices
        // virtual void resultantTangentMatrices
        // (
        //     // const vectorField& S,
        //     const diagTensorField& dSdE,
        //     tensor& DQDGamma,
        //     tensor& DQDK,
        //     tensor& DMDGamma,
        //     tensor& DMDK
        // ) const = 0;

        // //- Write cross-section scalar field in VTK format
        // virtual void writeVTK
        // (
        //     const fileName& fn,
        //     const scalarField& sf0,
        //     const scalarField& sf1,
        //     const vectorField& vf
        // ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
