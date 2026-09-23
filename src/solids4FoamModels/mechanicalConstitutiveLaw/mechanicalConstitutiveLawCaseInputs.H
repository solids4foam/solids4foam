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
    Foam::mechanicalConstitutiveLawCaseInputs

Description
    Scalar coupling inputs of mechanical constitutive laws that are read, at
    the current time, from another case directory.

    A law declares the scalar fields it reads through requiredScalarInputs(),
    and the manager normally takes each from the object registry, because
    another equation of the solid model solves for it. When nothing in this
    run solves for the field, a case can instead name a second case directory
    that holds it as a sequence of time directories, for example a
    temperature history computed beforehand by a thermal solver. That is what
    this class serves, for any input of any law.

    It is set per material, in the law's sub-dictionary of mechanicalProperties:

        steel
        {
            type    thermoMechanicalLaw;
            ...
            inputCaseDirectories
            {
                T   "hotCylinderTemperatureField";
            }
        }

    and the legacy spelling <input>caseDirectory, i.e. TcaseDirectory for the
    temperature of thermoMechanicalLaw, means the same, so an existing
    dictionary runs unchanged on either implementation. The path is relative
    to the case, as the legacy law's is.

    Semantics are those of the legacy thermoMechanicalLaw TcaseDirectory:
      - the second case has its own Time and mesh; the mesh must match this
        one (the same numbers of points, faces and cells and the same points),
        and additionally the same patches, since the boundary values are
        copied patch by patch;
      - the field is read once per time step, at the main case's current
        time, on the first evaluation of the step;
      - a step whose time directory holds no such field keeps the last one
        read, which is what lets a case give the field only at its start; a
        field never read at all is fatal;
      - the copy is registered under the input's name and written with the
        case's fields, for post-processing, unless something else already
        holds that name;
      - a field of the same name registered by another model, i.e. one this
        class did not read, takes precedence over the case directory.

    The copy lives on the full mesh. There are no per-material sub-meshes in
    the framework: each law gathers its own cells and boundary faces from it.
    Two materials may name different case directories for one input; each
    then has its own copy, and only the first is registered.

    In parallel the case directory is taken relative to this processor's case,
    exactly as the legacy law does: processorN/<caseDirectory> must hold that
    processor's portion of the second case, decomposed in the same way.

SourceFiles
    mechanicalConstitutiveLawCaseInputs.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawCaseInputs_H
#define mechanicalConstitutiveLawCaseInputs_H

#include "fvMesh.H"
#include "volFields.H"
#ifdef FOAMEXTEND
    #include "foamTime.H"
#else
    #include "Time.H"
#endif
#include "autoPtr.H"
#include "PtrList.H"
#include "HashTable.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class mechanicalConstitutiveLawCaseInputs Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawCaseInputs
{
    // Private Data

        //- One second case directory, with its own Time and mesh, both
        //  constructed on first use
        struct sourceCase
        {
            fileName caseDir_;
            autoPtr<Time> runTimePtr_;
            autoPtr<fvMesh> meshPtr_;
        };

        //- One input read from one case directory, copied onto the mesh
        struct sourcedField
        {
            //- Index into cases_
            label caseI_;

            //- Field name, which is the input's name
            word name_;

            //- The copy on the main mesh; empty until first read
            autoPtr<volScalarField> fieldPtr_;

            //- Time index the copy was last refreshed for
            label timeIndex_;
        };

        //- The main mesh
        const fvMesh& mesh_;

        //- Second cases, one per distinct directory
        mutable PtrList<sourceCase> cases_;

        //- Sourced fields, one per distinct (directory, name) pair
        mutable PtrList<sourcedField> fields_;

        //- Per law: input name -> index into fields_
        List<HashTable<label>> lawSources_;

        //- Inputs already reported as shadowed by a registered field
        mutable HashSet<word> reportedShadowed_;


    // Private Member Functions

        //- Construct the second case's Time and mesh and check the mesh
        void makeCase(sourceCase& sc) const;

        //- Check the second case's mesh against the main one. The full check,
        //  including the points, is done once; the sizes are checked on every
        //  read so that a topology change of the main mesh is refused rather
        //  than read out of bounds
        void checkMesh(const sourceCase& sc, const bool checkPoints) const;

        //- Read the field at the current time if not yet done this step
        void refresh(sourcedField& sf) const;

        //- No copy construct
        mechanicalConstitutiveLawCaseInputs
        (
            const mechanicalConstitutiveLawCaseInputs&
        ) = delete;

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawCaseInputs&) = delete;


public:

    // Static Member Functions

        //- Read the case-directory sources a law's dictionary gives for its
        //  inputs, in either spelling, refusing any that name an input the
        //  law does not read. Returns input name -> directory
        static HashTable<fileName> readSources
        (
            const dictionary& lawDict,
            const word& lawName,
            const wordList& requiredInputs
        );


    // Constructors

        //- Construct for a number of laws
        mechanicalConstitutiveLawCaseInputs
        (
            const fvMesh& mesh,
            const label nLaws
        );


    //- Destructor
    ~mechanicalConstitutiveLawCaseInputs();


    // Member Functions

        //- Source input 'name' of law lawI from case directory caseDir
        void addSource
        (
            const label lawI,
            const word& name,
            const fileName& caseDir
        );

        //- True if input 'name' of law lawI is read from a case directory
        bool found(const label lawI, const word& name) const;

        //- True if this object is one of the copies held here
        bool owns(const regIOobject& obj) const;

        //- The field for input 'name' of law lawI at the current time.
        //  Reads it from the case directory on the first call of a step.
        //  The manager is const where it gathers, so this is too, and the
        //  refresh is a cache update
        const volScalarField& field(const label lawI, const word& name) const;

        //- Say, once per input, that a registered field is used in place of
        //  the case directory
        void reportShadowed(const label lawI, const word& name) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
