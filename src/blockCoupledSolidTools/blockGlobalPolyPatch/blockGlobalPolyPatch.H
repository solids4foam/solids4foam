/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    blockGlobalPolyPatch

Description
    Global patch for block coupled solid solver.

Author
    Philip Cardiff UCD

SourceFiles
    blockGlobalPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef blockGlobalPolyPatch_H
#define blockGlobalPolyPatch_H

#include "coupledPolyPatch.H"
#include "transform.H"
#include "word.H"
#include "faceZone.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class blockGlobalPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class blockGlobalPolyPatch
:
    public coupledPolyPatch
{
protected:

    // Protected Member functions

        //- Initialise the calculation of the patch addressing
        virtual void initAddressing();

        //- Calculate the patch addressing
        virtual void calcAddressing();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();

public:

    //- Runtime type information
    TypeName("blockGlobalSolid");


    // Constructors

        //- Construct from components
        blockGlobalPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        blockGlobalPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        blockGlobalPolyPatch
        (
            const blockGlobalPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        blockGlobalPolyPatch
        (
            const blockGlobalPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new blockGlobalPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new blockGlobalPolyPatch
                (
                    refCast<const blockGlobalPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~blockGlobalPolyPatch();


    // Member functions

        // Patch ordering

            //- Initialize ordering for primitivePatch. Does not
            //  refer to *this (except for name() and type() etc.)
            virtual void initOrder(const primitivePatch&) const;

            //- Return new ordering for primitivePatch.
            //  Ordering is -faceMap: for every face
            //  index of the new face -rotation: for every new face the
            //  clockwise shift of the original face. Return false if nothing
            //  changes (faceMap is identity, rotation is 0), true otherwise.
            virtual bool order
            (
                const primitivePatch&,
                labelList& faceMap,
                labelList& rotation
            ) const;

            //- Synchronise communications of ordering for primitivePatch
            //  Used in cases when no topological change happens locally,
            //  but is happening on other processors
            virtual void syncOrder() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
