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
    ExtendedGGIInterpolation

Description
    Mass-conservative face interpolation of face data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

SourceFiles
    ExtendedGGIInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef ExtendedGGIInterpolation_H
#define ExtendedGGIInterpolation_H

//#include "GGIInterpolationTemplate.H"
#include "newGGIInterpolationTemplate.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ExtendedGGIInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class MasterPatch, class SlavePatch>
class ExtendedGGIInterpolation
:
    //public GGIInterpolation<MasterPatch, SlavePatch>
    public newGGIInterpolation<MasterPatch, SlavePatch>
{
    // Private data


    // Demand-driven data

        //- Slave face into which each point of master patch is projected
        mutable List<labelPair>* masterPointAddressingPtr_;

        //- Weighting factors
        mutable FieldField<Field, scalar>* masterPointWeightsPtr_;

        //- Distance to intersection for master patch points
        mutable scalarField* masterPointDistancePtr_;

        //- Master face into which each point of slave patch is projected
        mutable List<labelPair>* slavePointAddressingPtr_;

        //- Weighting factors
        mutable FieldField<Field, scalar>* slavePointWeightsPtr_;

        //- Distance to intersection for master patch points
        mutable scalarField* slavePointDistancePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ExtendedGGIInterpolation(const ExtendedGGIInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const ExtendedGGIInterpolation&);

        //- Calculate point weights
        void calcMasterPointAddressing() const;

        //- Calculate point weights
        void calcMasterPointWeights() const;

        //- Calculate point weights
        void calcSlavePointAddressing() const;

        //- Calculate point weights
        void calcSlavePointWeights() const;

    public:

    // Static data

        ClassName("ExtendedGGIInterpolation");

    // Constructors

        //- Construct from components
        ExtendedGGIInterpolation
        (
            const MasterPatch& masterPatch,
            const SlavePatch&  slavePatch,
            const tensorField& forwardT,
            const tensorField& reverseT,
            const vectorField& forwardSep,
            const scalar masterFaceNonOverlapFaceTol = 0,
            const scalar slaveFaceNonOverlapFaceTol = 0,
            const bool rescaleGGIWeightingFactors = true,
            // const GGIInterpolationName::quickReject reject =
            //     GGIInterpolationName::AABB,
            const newGGIInterpolationName::quickReject reject =
                newGGIInterpolationName::AABB //,
            // const boundBox& regionOfInterest =
            //     boundBox(vector::min, vector::max)
        );


    // Destructor

        ~ExtendedGGIInterpolation();


    // Member Functions

        // Access

        //- Return reference to point addressing
        const List<labelPair>& masterPointAddr() const;

        //- Return reference to point weights
        const FieldField<Field, scalar>& masterPointWeights() const;

        //- Return distance to intersection for patch points
        const scalarField& masterPointDistanceToIntersection() const;

        //- Return reference to point addressing
        const List<labelPair>& slavePointAddr() const;

        //- Return reference to point weights
        const FieldField<Field, scalar>& slavePointWeights() const;

        //- Return distance to intersection for patch points
        const scalarField& slavePointDistanceToIntersection() const;

    // Interpolation functions

        //- Interpolate point field
        template<class Type>
        tmp<Field<Type> > masterToSlavePointInterpolate
        (
            const Field<Type>& pf
        ) const;

        //         template<class Type>
        //         tmp<Field<Type> > masterToSlavePointInterpolate
        //         (
        //             const tmp<Field<Type> >& tpf
        //         ) const;

        template<class Type>
        tmp<Field<Type> > slaveToMasterPointInterpolate
        (
            const Field<Type>& pf
        ) const;

        //         template<class Type>
        //         tmp<Field<Type> > slaveToMasterPointInterpolate
        //         (
        //             const tmp<Field<Type> >& tpf
        //         ) const;

    // Edit

        //- Correct weighting factors for moving mesh.
        //bool movePoints();
        bool movePoints
        (
            const tensorField& forwardT,
            const tensorField& reverseT,
            const vectorField& forwardSep
        );
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "ExtendedGGIInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
