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
    Foam::mechanicalConstitutiveLawStateIO

Description
    Reading and writing the persistent state of a mechanical constitutive law,
    so that a history dependent material can restart.

    Two things want a state field written, and they want different things:

      - a restart needs the value back exactly as it was, otherwise it
        continues a different calculation that merely resembles the first;

      - a user wants to look at it, and ParaView can show a cell field or a
        point field and nothing else.

    Trying to serve both with one artefact is what makes this awkward. They are
    separated here. This class does the first: it writes the state in the shape
    the state is already held, a flat Field, as an IOField. That works for every
    integration point topology, including the quadrature ones where the points
    correspond to no mesh entity at all and no GeometricField could be formed.
    The second is a projection, written and never read, and does not belong
    here.

    The cost of an IOField is that decomposePar and reconstructPar do not map
    it, so the state can only be read back on the decomposition that wrote it.
    That is a real limit, and the alternative was worse: a volField restart
    would buy decomposition changes for two of the six topologies at the price
    of two IO paths, a dimensionSet that no law declares, and a boundary field
    whose evaluation would overwrite the very history being restored. A limit
    that announces itself is preferable to a silent approximation, so a
    mismatch here is fatal rather than mapped.

    One file per law, per topology, per state variable, holding that law's
    integration points followed by its boundary points patch by patch. The name
    carries the law and the topology because two laws may legitimately hold
    different values of the same variable name, and the same law may hold state
    on more than one topology in the same run.

SourceFiles
    mechanicalConstitutiveLawStateIO.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawStateIO_H
#define mechanicalConstitutiveLawStateIO_H

#include "regIOobject.H"
#include "mechanicalConstitutiveLawState.H"
#include "mechanicalConstitutiveLawStateSpec.H"
#include "fvMesh.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class mechanicalConstitutiveLawStateIO Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawStateIO
{
public:

    //- The pieces of one state variable, in the order they are written: the
    //  law's own integration points, then its points on each patch. Held as
    //  states rather than as fields because a state may reallocate its fields
    //  and the state objects themselves do not move
    typedef List<mechanicalConstitutiveLawState*> stateParts;


    // Static Member Functions

        //- Name of the file holding one state variable
        //  law, topology and the path through any child states, so that two
        //  laws declaring the same variable name cannot collide
        static word fieldName
        (
            const word& lawName,
            const word& topologyName,
            const wordList& childPath,
            const word& variableName
        );

        //- Total number of entries across the parts
        static label totalSize(const stateParts& parts);

        //- Which decomposition this state belongs to, as a short string
        //  carried in the file's header note.
        //
        //  A flat list means something only against the mesh it was written
        //  for. Two decompositions of the same case into the same number of
        //  processors can give a rank the same *number* of cells and a
        //  different set of them, and then a stale processor directory left
        //  behind by an earlier decomposition holds a file of exactly the
        //  right length and entirely the wrong values. A length check cannot
        //  see that; this can
        static string decompositionIdentity(const fvMesh& mesh);

        //- The cell each of this processor's cells was in the serial mesh,
        //  or an empty list when the case was not decomposed
        static labelList cellProcAddressing(const fvMesh& mesh);

        //- The face each of this processor's faces was in the serial mesh.
        //  decomposePar stores these offset by one and signed, so that face
        //  zero can carry a sign and a reversed face can be told from a
        //  forward one; only the magnitude matters here, since state belongs
        //  to a face rather than to a direction across it
        static labelList faceProcAddressing(const fvMesh& mesh);

        //- Path of the undecomposed case, which is where the serial state
        //  sits when a decomposed run has to go looking for it
        static fileName globalCasePath(const fvMesh& mesh);

        //- The serial cell count recorded in a serial file's header note
        static label serialCellCount(const string& note);

        //- How an integration point says which piece of mesh it belongs to.
        //  Cells and faces share one numbering so that a single list can carry
        //  a law's own points and its boundary points together: a cell is
        //  itself, a face is offset past every cell
        static inline label cellEntity(const label cellI)
        {
            return cellI;
        }

        static inline label faceEntity(const label nCells, const label faceI)
        {
            return nCells + faceI;
        }
};


/*---------------------------------------------------------------------------*\
                      Struct stateFieldAccess Declaration
\*---------------------------------------------------------------------------*/

//- Reaching a state field of a given type by name.
//  The state exposes one accessor per type rather than a template, so this
//  turns the four of them back into one, and the IO below is written once
template<class Type>
struct stateFieldAccess;

template<>
struct stateFieldAccess<scalar>
{
    static const Field<scalar>& get
    (
        const mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.getScalarField(n);
    }

    static Field<scalar>& ref
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.scalarField(n);
    }

    static Field<scalar>& ref0
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.scalarField0(n);
    }
};

template<>
struct stateFieldAccess<vector>
{
    static const Field<vector>& get
    (
        const mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.getVectorField(n);
    }

    static Field<vector>& ref
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.vectorField(n);
    }

    static Field<vector>& ref0
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.vectorField0(n);
    }
};

template<>
struct stateFieldAccess<tensor>
{
    static const Field<tensor>& get
    (
        const mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.getTensorField(n);
    }

    static Field<tensor>& ref
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.tensorField(n);
    }

    static Field<tensor>& ref0
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.tensorField0(n);
    }
};

template<>
struct stateFieldAccess<symmTensor>
{
    static const Field<symmTensor>& get
    (
        const mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.getSymmTensorField(n);
    }

    static Field<symmTensor>& ref
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.symmTensorField(n);
    }

    static Field<symmTensor>& ref0
    (
        mechanicalConstitutiveLawState& s, const word& n
    )
    {
        return s.symmTensorField0(n);
    }
};


/*---------------------------------------------------------------------------*\
                    Class stateIOFieldProxy Declaration
\*---------------------------------------------------------------------------*/

//- Writes one state variable when the registry writes.
//  The state itself stays the authoritative copy: this gathers from it at
//  write time rather than holding a duplicate that could drift. writeData is
//  overridden rather than writeObject because writeObject's signature differs
//  between the three forks and writeData's does not
template<class Type>
class stateIOFieldProxy
:
    public regIOobject
{
    // Private data

        //- The state objects to gather from, in write order
        const mechanicalConstitutiveLawStateIO::stateParts parts_;

        //- Name of the state variable within each state
        const word variableName_;


public:

    // Constructors

        stateIOFieldProxy
        (
            const IOobject& io,
            const mechanicalConstitutiveLawStateIO::stateParts& parts,
            const word& variableName,
            const fvMesh& mesh
        )
        :
            regIOobject(io),
            parts_(parts),
            variableName_(variableName)
        {
            // Stamped into the header so that a restart can tell this file
            // apart from one of the same length written for a different
            // decomposition
            note() =
                mechanicalConstitutiveLawStateIO::decompositionIdentity(mesh);
        }


    //- Destructor
    virtual ~stateIOFieldProxy() = default;


    // Member Functions

        //- Announce the same class an IOField would.
        //  The header carries this, and the read side is a plain IOField: a
        //  proxy that named itself would write a file only it could read
        virtual const word& type() const
        {
            return IOField<Type>::typeName;
        }

        //- Gather the state and write it as a plain list
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "mechanicalConstitutiveLawStateIOTemplates.C"
#endif

#endif

// ************************************************************************* //
