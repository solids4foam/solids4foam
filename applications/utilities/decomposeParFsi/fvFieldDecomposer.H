/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fvFieldDecomposer

Description
    Finite Volume volume and surface field decomposer.

SourceFiles
    fvFieldDecomposer.C
    fvFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef fvFieldDecomposer_H
#define fvFieldDecomposer_H

#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class fvFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class fvFieldDecomposer
{
public:

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            // Private data

                label sizeBeforeMapping_;
                labelList directAddressing_;

        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const label sizeBeforeMapping,
                    const unallocLabelList& addressingSlice,
                    const label addressingOffset
                );


            // Member functions

                label size() const
                {
                    return directAddressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return true;
                }

                const unallocLabelList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class
        class processorVolPatchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            // Private data

                label sizeBeforeMapping_;
                labelListList addressing_;
                scalarListList weights_;

        public:

            //- Construct given addressing
            processorVolPatchFieldDecomposer
            (
                const fvMesh& mesh,
                const unallocLabelList& addressingSlice
            );


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


        //- Processor patch field decomposer class
        class processorSurfacePatchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            label sizeBeforeMapping_;
            labelListList addressing_;
            scalarListList weights_;

        public:

            //- Construct given addressing
            processorSurfacePatchFieldDecomposer
            (
                label sizeBeforeMapping,
                const unallocLabelList& addressingSlice
            );


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


private:

    // Private data

        //- Reference to complete mesh
        const fvMesh& completeMesh_;

        //- Reference to processor mesh
        const fvMesh& procMesh_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        List<patchFieldDecomposer*> patchFieldDecomposerPtrs_;

        List<processorVolPatchFieldDecomposer*>
            processorVolPatchFieldDecomposerPtrs_;

        List<processorSurfacePatchFieldDecomposer*>
            processorSurfacePatchFieldDecomposerPtrs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvFieldDecomposer(const fvFieldDecomposer&);

        //- Disallow default bitwise assignment
        void operator=(const fvFieldDecomposer&);


public:

    // Constructors

        //- Construct from components
        fvFieldDecomposer
        (
            const fvMesh& completeMesh,
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );


    // Destructor

        ~fvFieldDecomposer();


    // Member Functions

        //- Decompose volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        decomposeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        decomposeField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
