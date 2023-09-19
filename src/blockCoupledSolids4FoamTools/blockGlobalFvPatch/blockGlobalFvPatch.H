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
    blockGlobalFvPatch

Description
    Global patch for block coupled solid solver.

Author
    Philip Cardiff UCD

SourceFiles
    blockGlobalFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef blockGlobalFvPatch_H
#define blockGlobalFvPatch_H

#include "coupledFvPatch.H"
#include "blockGlobalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class blockGlobalFvPatch Declaration
\*---------------------------------------------------------------------------*/

class blockGlobalFvPatch
:
    public coupledFvPatch
    //public ggiLduInterface
{
    // Private Data

        //- Reference to polyPatch
        const blockGlobalPolyPatch& blockGlobalPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(vectorField&) const;


public:

    //- Runtime type information
    TypeName(blockGlobalPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        blockGlobalFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            blockGlobalPolyPatch_(refCast<const blockGlobalPolyPatch>(patch))
        {}


    // Destructor

        virtual ~blockGlobalFvPatch();


    // Member functions

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
