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
    Foam::mechanicalConstitutiveLawInputGatherer

Description
    Gathers the scalar coupling inputs each mechanical constitutive law of a
    mechanicalConstitutiveLawManager reads, such as a temperature, onto the
    law's integration points.

    An input is taken from the object registry, where another model solves
    for it; failing that, from another case directory, as
    mechanicalConstitutiveLawCaseInputs serves it; and failing that, from a
    file in the current or initial time directory, read once per time step on
    every rank. It also warns, once per field, where a law could read a
    registered field but was not asked to.

    Split out of the manager, which holds one and forwards to it. It holds
    references to the manager's mesh, laws, law cells, case inputs and
    reported-input set, and owns the cache of inputs read from files.

SourceFiles
    mechanicalConstitutiveLawInputGatherer.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawInputGatherer_H
#define mechanicalConstitutiveLawInputGatherer_H

#include "fvMesh.H"
#include "volFields.H"
#include "HashSet.H"
#include "mechanicalConstitutiveLaw.H"
#include "integrationPointTopology.H"
#include "mechanicalConstitutiveLawCaseInputs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class mechanicalConstitutiveLawInputGatherer Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawInputGatherer
{
    // Private Data

        //- The mesh the laws are defined on
        const fvMesh& mesh_;

        //- The laws, as the manager holds them
        const PtrList<mechanicalConstitutiveLaw>& laws_;

        //- The cells of each law
        const List<labelList>& lawCells_;

        //- Scalar inputs read from another case directory; empty if no law
        //  asks for that
        const autoPtr<mechanicalConstitutiveLawCaseInputs>& caseInputsPtr_;

        //- The optional coupling inputs already reported as registered but
        //  unread. Shared with the manager, which reports unread prescribed
        //  fields through the same set
        HashSet<word>& reportedUnreadInputs_;

        //- Coupling inputs that exist only as files, read once per time step
        //  on every rank by refreshScalarInputs, and the time index they were
        //  read at. Unregistered, so that nothing else takes them for a field
        //  another model supplies
        mutable HashTable<autoPtr<volScalarField>> diskScalarInputs_;
        mutable label diskScalarInputsTimeIndex_;


    // Private Member Functions

        //- The field a law's scalar input is gathered from, or nullptr if
        //  neither the registry nor a case directory supplies it, in which
        //  case the caller falls back to a field on disk.
        //
        //  A field registered by another model comes first, so that one
        //  solved for in this run is always what the law sees. A case
        //  directory source comes next; the copy it registers is recognised
        //  as its own and refreshed, once per time step, rather than taken as
        //  another model's
        const volScalarField* scalarInputSource
        (
            const label lawI,
            const word& name
        ) const;

        //- A coupling input that exists only as a file, as refreshScalarInputs
        //  read it for this time step, or nullptr if there is none
        const volScalarField* diskScalarInput(const word& name) const;


public:

    // Constructors

        //- Construct from references to the manager's members, which must
        //  outlive this
        mechanicalConstitutiveLawInputGatherer
        (
            const fvMesh& mesh,
            const PtrList<mechanicalConstitutiveLaw>& laws,
            const List<labelList>& lawCells,
            const autoPtr<mechanicalConstitutiveLawCaseInputs>& caseInputsPtr,
            HashSet<word>& reportedUnreadInputs
        );


    // Member Functions

        //- Gather the live coupling inputs a law asked for, as views over its
        //  own integration points, and hand back an inputs object carrying
        //  them. A law that asks for nothing gets an inputs object holding
        //  only the time step, which is what every law had before
        mechanicalConstitutiveLawInputs lawInputs
        (
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            const scalar dt,
            List<HashTable<autoPtr<scalarField>>>& lawScalarInputs
        ) const;

        //- As lawInputs, for a law's faces on one boundary patch. The values
        //  come from the field's patch values, not from the cells behind it,
        //  and on a coupled patch are interpolated as on an internal face
        mechanicalConstitutiveLawInputs lawInputsPatch
        (
            const label lawI,
            const label patchI,
            const labelList& faces,
            const scalar dt,
            List<HashTable<autoPtr<scalarField>>>& lawScalarInputs
        ) const;

        //- Bring every coupling input source up to date for this time step:
        //  the case-directory inputs, and the inputs that exist only as files.
        //  Called at the start of every evaluation path, on every rank and
        //  before any material is skipped for having no points there, since
        //  reading a file can be collective
        void refreshScalarInputs() const;

        //- The coupling inputs this law needs, including every law below it.
        //
        //  A composite's sub-law reads its own inputs, and asking only the top
        //  of the tree loses them. State already recurses this way; this makes
        //  dependencies match, so that nesting one composite inside another
        //  gathers what the inner law asked for instead of silently omitting
        //  it
        static wordList requiredScalarInputsRecursive
        (
            const mechanicalConstitutiveLaw& law
        );

        //- Warn, once per field, where a law could read a registered field as
        //  a coupling input but was not asked to, and so uses its own
        //  substitute: the case then runs, but not on the field it has. This
        //  law and every law below it
        void reportUnreadScalarInputs
        (
            const mechanicalConstitutiveLaw& law
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
