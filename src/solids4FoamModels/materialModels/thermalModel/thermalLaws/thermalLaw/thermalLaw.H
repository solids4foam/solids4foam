/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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
    thermalLaw

Description
    Thermal material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    thermalLaw.C
    newThermalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef thermalLaw_H
#define thermalLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalLaw Declaration
\*---------------------------------------------------------------------------*/

class thermalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        thermalLaw(const thermalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const thermalLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

public:

    //- Runtime type information
    TypeName("thermalLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalLaw,
            dictionary,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<thermalLaw> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        thermalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~thermalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return specific heat capacity
        virtual tmp<volScalarField> C() const = 0;

        //- Return thermal conductivity
        virtual tmp<volScalarField> k() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
