/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    coulombFriction

Description
    coulombFriction friction law
    shear traction is linearly related to the contact pressure
    via the friction coefficient

SourceFiles
    coulombFriction.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef coulombFriction_H
#define coulombFriction_H

#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coulombFriction Declaration
\*---------------------------------------------------------------------------*/

class coulombFriction
:
    public frictionLaw
{
    // Private data

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- Coulomb friction coefficient
        scalar frictionCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        //coulombFriction(const coulombFriction&);

        //- Disallow default bitwise assignment
        void operator=(const coulombFriction&);

public:

    //- Runtime type information
    TypeName("coulomb");

    // Static data members


    // Constructors

        //- Construct from dictionary
        coulombFriction
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        coulombFriction(const coulombFriction& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new coulombFriction(*this)
                );
        }



    // Destructor

        virtual ~coulombFriction();


    // Member Functions

        //- Return slip traction field
        // the traction that would cause slipping for the given pressure
        virtual scalar slipTraction(const scalar pressure);

        //- Return slip traction field given the pressure and slip direction
        //  This allows for anisotropic friction
        virtual scalar slipTraction
        (
            const scalar pressure, const vector& slipDir
        );

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
