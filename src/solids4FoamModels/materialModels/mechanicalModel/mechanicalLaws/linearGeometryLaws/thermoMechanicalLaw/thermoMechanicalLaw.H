/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    thermoMechanicalLaw

Description
    Generalised thermo-coupled mechanical law, where a thermal volumetric stress
    term is add to the stress calculated from any other run-time selectable
    mechanical law.

    Careful: if the run-time selected mechanical law is a non-linear function of
    the volumetric stress then this de-coupled approach may not make sense. For
    example, if plasticity depended on the volumetric stress.

SourceFiles
    thermoMechanicalLaw.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Inspired by Denis Maier, BAW, code design suggestions.

\*---------------------------------------------------------------------------*/

#ifndef thermoMechanicalLaw_H
#define thermoMechanicalLaw_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class thermoMechanicalLaw
:
    public mechanicalLaw
{
    // Private Data

        //- Run-time selectable mechanical law to define the stress without a
        //  a thermal term
        autoPtr<mechanicalLaw> mechLawPtr_;

        //- Coefficient of linear expansion
        const dimensionedScalar alpha_;

        //- Stress-free reference temperature
        const dimensionedScalar T0_;

        //- If the temperature field T is not found in the objectRegistry then
        //  it is read from disk. This means that if the solidModel does not
        //  solve a temperature equation, then the user must specify a constant
        //  temperature field T
        mutable autoPtr<volScalarField> TPtr_;

        //- Flag to indicate the T field has be read from disk at least once
        bool TFieldWasReadFromDisk_;

        //- Relative address of the temperature case
        const fileName TcaseDir_;

        //- Time object for the temperature case
        mutable autoPtr<Time> TrunTimePtr_;

        //- fvMesh object for the temperature case
        mutable autoPtr<fvMesh> TmeshPtr_;

        //- Store the current time index so we know the start of a new time-step
        label curTimeIndex_;

    // Private Member Functions

        //- Make the TrunTime object
        void makeTrunTime();

        //- Return a reference to TrunTime
        Time& TrunTime();

        //- Make the Tmesh object
        void makeTmesh();

        //- Return a reference to Tmesh
        fvMesh& Tmesh();

        //- Read the temperature field from disk
        //  Return true if the T field has ever been read from disk
        bool readTField();

        //- Return a reference to the initial pore-temperature surface field
        const surfaceScalarField& p0f() const;

        //- Lookup and return a reference to the temperature field
        const volScalarField& lookupTemperatureField();

        //- Disallow default bitwise copy construct
        thermoMechanicalLaw(const thermoMechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const thermoMechanicalLaw&);

public:

    //- Runtime type information
    TypeName("thermoMechanicalLaw");

    // Static data members


    // Constructors

        //- Construct from dictionary
        thermoMechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~thermoMechanicalLaw();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

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
