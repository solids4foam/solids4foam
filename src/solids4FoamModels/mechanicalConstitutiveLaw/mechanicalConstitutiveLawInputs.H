/*---------------------------------------------------------------------------*\
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

Class
    Foam::mechanicalConstitutiveLawInputs

Description
    Immutable, per-evaluation inputs handed to a mechanical constitutive law
    alongside its kinematics, state and response.

    These are values a law reads and never writes, which are neither
    kinematics nor history: the time increment, and fields owned by the solid
    model and solved by some other equation, such as a temperature, a pore
    pressure, or the pressure of a mixed displacement-pressure formulation.

    Why this is a separate object from the kinematics, rather than more
    members on it:

      - The kinematics type is the formulation boundary. There are two
        evaluate() overloads, and a finite-strain-only law refuses a
        small-strain evaluation by not implementing that one. A temperature
        has no business taking part in that choice.

      - The finite-difference tangent reconstructs the kinematics object once
        per Voigt component, forwarding by hand everything it is not
        perturbing. Anything living in the kinematics is therefore one more
        thing to forget to forward, silently. This object is passed through
        each perturbation unchanged, so it cannot be dropped.

    Why it is not carried in the state, which was the first design:

        A shadow state aliases the parent's old-time fields and owns its own,
        fresh, current-time fields, which is what makes a tangent query safe.
        A live input written into current-time state would therefore be absent
        inside the shadow, so every tangent query and every finite-difference
        perturbation would evaluate with a missing input. Putting it in
        old-time state is worse: it disguises a current input as history, and
        storeOldTime() would commit it.

    Fields are supplied per integration point, indexed exactly as the law's
    own evaluation loop is, so a law reads inputs.getScalar("T")[i] beside its
    kinematics. The law never sees a mesh, a registry or a GeometricField.

    Absence is not zero. A law that can work without an input asks with
    findScalar(), which returns nullptr; one that requires an input asks with
    getScalar(), which fails rather than silently reading zero.

SourceFiles
    mechanicalConstitutiveLawInputs.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawInputs_H
#define mechanicalConstitutiveLawInputs_H

#include "scalarField.H"
#include "vectorField.H"
#include "HashTable.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mechanicalConstitutiveLawInputs Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawInputs
{
    // Private data

        //- Time increment
        const scalar dt_;

        //- Scalar fields, one value per integration point of the law being
        //  evaluated. Not owned: the supplier outlives the evaluation
        HashTable<const scalarField*> scalarFields_;

        //- Vector fields, likewise
        HashTable<const vectorField*> vectorFields_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawInputs&) = delete;


public:

    // Constructors

        //- Construct from the time increment
        explicit mechanicalConstitutiveLawInputs(const scalar dt)
        :
            dt_(dt),
            scalarFields_(),
            vectorFields_()
        {}


    //- Destructor
    ~mechanicalConstitutiveLawInputs() = default;


    // Member Functions

        // Access

            //- Time increment
            scalar dt() const
            {
                return dt_;
            }

            //- Is a scalar input available under this name?
            bool foundScalar(const word& name) const
            {
                return scalarFields_.found(name);
            }

            //- Is a vector input available under this name?
            bool foundVector(const word& name) const
            {
                return vectorFields_.found(name);
            }

            //- Scalar input, or nullptr where it was not supplied.
            //  For a law that can work without it
            const scalarField* findScalar(const word& name) const;

            //- Vector input, or nullptr where it was not supplied
            const vectorField* findVector(const word& name) const;

            //- Scalar input, failing where it was not supplied.
            //  For a law that requires it: a missing input must not read as
            //  zero, which would be a plausible wrong answer
            const scalarField& getScalar(const word& name) const;

            //- Vector input, failing where it was not supplied
            const vectorField& getVector(const word& name) const;


        // Edit

            //- Supply a scalar input. The field must outlive this object and
            //  must have one value per integration point of the law
            void setScalar(const word& name, const scalarField& fld);

            //- Supply a vector input
            void setVector(const word& name, const vectorField& fld);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
