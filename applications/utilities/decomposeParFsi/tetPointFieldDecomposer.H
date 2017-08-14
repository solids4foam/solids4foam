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
    tetPointFieldDecomposer

Description
    Tetrahedral point field decomposer.

SourceFiles
    tetPointFieldDecomposer.C
    tetPointFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef tetPointFieldDecomposer_H
#define tetPointFieldDecomposer_H

#include "tetPolyMesh.H"
#include "tetPointMesh.H"
#include "elementMesh.H"
#include "tetPolyPatch.H"
#include "PointPatchFieldMapperPatchRef.H"
#include "tetPolyPatchFields.H"
#include "elementPatchFields.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class tetPointFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class tetPointFieldDecomposer
{
public:

        //- Patch field decomposer class
        class tetPolyPatchFieldDecomposer
        :
            public PointPatchFieldMapperPatchRef<tetPolyPatch>
        {
            // Private data

                //- Size before mapping
                label sizeBeforeMapping_;

                //- Reference to direct addressing
                const labelList& directAddressing_;


            // Demand-driven data

                //- Addressing
                mutable labelList* directPatchAddressingPtr_;


            // Private Member Functions

                //- Disallow default bitwise copy construct
                tetPolyPatchFieldDecomposer(const tetPolyPatchFieldDecomposer&);

                //- Disallow default bitwise assignment
                void operator=(const tetPolyPatchFieldDecomposer&);


            // Member functions to calculate demand driven data

                //- Calculate addressing
                void calcPatchAddressing() const;


        public:

            // Constructors

                //- Construct given addressing
                tetPolyPatchFieldDecomposer
                (
                    const tetPolyPatch& sourcePatch,
                    const tetPolyPatch& targetPatch,
                    const labelList& directAddr
                )
                :
                    PointPatchFieldMapperPatchRef<tetPolyPatch>
                    (
                        sourcePatch,
                        targetPatch
                    ),
                    sizeBeforeMapping_(sourcePatch.size()),
                    directAddressing_(directAddr),
                    directPatchAddressingPtr_(NULL)
                {}


            // Destructor

                virtual ~tetPolyPatchFieldDecomposer()
                {}


            // Member functions

                virtual label size() const
                {
                    return directAddressing().size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                virtual bool direct() const
                {
                    return true;
                }

                virtual const unallocLabelList& directAddressing() const
                {
                    if (!directPatchAddressingPtr_)
                    {
                        calcPatchAddressing();
                    }

                    return *directPatchAddressingPtr_;
                }
        };

private:

    // Private data

        //- Reference to original mesh
        const tetPolyMesh& originalMesh_;

        //- Reference to processor mesh
        const tetPolyMesh& processorMesh_;

        //- Reference to point addressing
        const labelList& pointAddressing_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        List<tetPolyPatchFieldDecomposer*> patchFieldDecompPtrs_;

    // Demand-driven data

        //- Addressing
        mutable labelList* directAddressingPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPointFieldDecomposer
        (
            const tetPointFieldDecomposer&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPointFieldDecomposer&);

    // Member functions for demand driven data

        //- Return addressing
        const labelList& directAddressing() const;

        //- Calculate addressing
        void calcAddressing() const;


public:

    // Constructors

        //- Construct from components
        tetPointFieldDecomposer
        (
            const tetPolyMesh& originalMesh,
            const tetPolyMesh& processorMesh,
            const labelList& pointAddressing,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );

    // Destructor

        ~tetPointFieldDecomposer();


    // Member Functions

        //- Decompose point field
        template<class Type>
        tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
        decomposeField
        (
            const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        ) const;

        //- Decompose element field
        template<class Type>
        tmp<GeometricField<Type, elementPatchField, elementMesh> >
        decomposeField
        (
            const GeometricField<Type, elementPatchField, elementMesh>&
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetPointFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
