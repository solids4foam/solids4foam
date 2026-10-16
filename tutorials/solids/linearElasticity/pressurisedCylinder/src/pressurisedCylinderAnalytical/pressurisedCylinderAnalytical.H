/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    solids4foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    pressurisedCylinderAnalytical

Description
    Analytical displacement and stress for the internally pressurised cylinder.

    Analytical solution taken from:
    J. G. S.P. Timoshenko, Theory of Elasticity. McGraw-Hill, New York, 1970.

SourceFiles
    pressurisedCylinderAnalytical.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pressurisedCylinderAnalytical_H
#define pressurisedCylinderAnalytical_H

#include "dictionary.H"
#include "point.H"
#include "vector.H"
#include "vectorField.H"
#include "symmTensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                       Class pressurisedCylinderAnalytical Declaration
\*---------------------------------------------------------------------------*/

class pressurisedCylinderAnalytical
{
    // Private data

        //- Inside cylinder radius
        scalar Ri_;

        //- Outside cylinder radius
        scalar Ro_;

        //- Applied internal pressure (external pressure is zero)
        scalar p_;

        //- Young's modulus
        scalar E_;

        //- Poisson's ratio
        scalar nu_;

        //- Helper value (p_*sqr(Ri_)/(sqr(Ro_) - sqr(Ri_)))
        scalar A_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const pressurisedCylinderAnalytical&);

public:

    //- Runtime type information
    TypeName("pressurisedCylinderAnalytical");


    // Constructors

        //- Empty constructor
        pressurisedCylinderAnalytical();

        //- Copy constructor
        pressurisedCylinderAnalytical(const pressurisedCylinderAnalytical& pc);

        //- Construct from a dictionary
        pressurisedCylinderAnalytical(const dictionary& dict);


    // Destructor

        virtual ~pressurisedCylinderAnalytical()
        {}


    // Member Functions

        // Access

            //- Return the displacements at the given location
            vector displacement(const vector& location) const;

            //- Return the displacements (in cylindrical coordinate system)
            //  at the given location
            vector cylindricalDisplacement(const vector& location) const;

            //- Return the displacements at the given locations
            tmp<vectorField> displacement(const vectorField& locations) const;

            //- Return the stress tensor at the given location
            symmTensor stress(const vector& location) const;

            //- Return the stress tensor (in cylindrical coordinate system)
            //  at the given location
            symmTensor cylindricalStress(const vector& location) const;

    // Write

        //- Write parameter dict to the given stream
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
