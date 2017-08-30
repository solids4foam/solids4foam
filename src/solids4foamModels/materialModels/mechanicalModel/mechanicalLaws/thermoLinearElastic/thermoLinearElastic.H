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
    thermoLinearElastic

Description
    Duhamel-Neumann form of Hooke's law for linear thermo-elasticity.

SourceFiles
    thermoLinearElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef thermoLinearElastic_H
#define thermoLinearElastic_H

#include "linearElastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoLinearElastic Declaration
\*---------------------------------------------------------------------------*/

class thermoLinearElastic
:
    public linearElastic
{
    // Private data

        //- Coefficient of linear expansion
        const dimensionedScalar alpha_;

        //- Stress-free reference temperature
        const dimensionedScalar T0_;

        //- If the temperature field T is not found in the objectRegistry then
        //  it is read from disk. This means that if the solidModel does not
        //  solve a temperature equation, then the user must specify a constant
        //  temperature field T
        mutable autoPtr<volScalarField> TPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoLinearElastic(const thermoLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const thermoLinearElastic&);

public:

    //- Runtime type information
    TypeName("thermoLinearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        thermoLinearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~thermoLinearElastic();


    // Member Functions

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
