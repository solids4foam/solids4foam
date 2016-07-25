/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    mechanicalModel

Description
    Class in charge of mechanical properties

SourceFiles
    mechanicalModel.C

Author
    Philip Cardiff UCD
    philip.cardiff@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef mechanicalModel_H
#define mechanicalModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "mechanicalLaw.H"
#include "Switch.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class mechanicalModel Declaration
\*---------------------------------------------------------------------------*/

class mechanicalModel
:
  public IOdictionary
{
    // Private data

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- mechanical law
        autoPtr<mechanicalLaw> mechanicalLawPtr_;

        //- Plane stress
        Switch planeStress_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mechanicalModel(const mechanicalModel&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalModel&);


public:

    //- Runtime type information
    TypeName("mechanicalModel");


    // Constructors

        //- Construct
        mechanicalModel(const fvMesh& mesh);


    // Destructor

        virtual ~mechanicalModel();


    // Member Functions

        // Access

        //- Return const reference to mesh
        const fvMesh& mesh() const;

        //- Return true for plane stress
        const Switch& planeStress() const;

        //- Return mechanical law
        const mechanicalLaw& law() const;

        //- Return density
        tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        tmp<volScalarField> impK() const;

        //- Return the implicit stiffness surface field
        //  This is the diffusivity for the Laplacian term
        tmp<surfaceScalarField> impKf() const;

        //- Update stress based on plasticity stress return vol field
        virtual void correct(volSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update yield stress
        void updateYieldStress();

        //- Rate of dissipated energy from plasticity
        //tmp<volScalarField> plasticDissipationRate() const;

        //- Read plasticityProperties dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
