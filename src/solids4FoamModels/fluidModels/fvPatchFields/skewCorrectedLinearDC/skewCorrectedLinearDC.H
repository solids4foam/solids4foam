/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::skewCorrectedLinearDC

Description
    skewCorrectedLinearDC interpolation scheme uses extrapolated gradient-based
    correction from both sides onto the face, using 1/2 weighting.
    Class is derived from upwind and returns upwind weighting factors.
    Correction is calculated using difference between linear and upwind interp.

Author
    Zeljko Tukovic, FSB Zagreb  All rights reserved.

SourceFiles
    skewCorrectedLinearDC.C

\*---------------------------------------------------------------------------*/

#ifndef skewCorrectedLinearDC_H
#define skewCorrectedLinearDC_H

#include "upwind.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class skewCorrectedLinearDC Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewCorrectedLinearDC
:
    public upwind<Type>
{
    // Private Data

        tmp<fv::gradScheme<Type> > gradScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewCorrectedLinearDC(const skewCorrectedLinearDC&);

        //- Disallow default bitwise assignment
        void operator=(const skewCorrectedLinearDC&);


public:

    //- Runtime type information
    TypeName("skewCorrectedLinearDC");


    // Constructors

        //- Construct from faceFlux
        skewCorrectedLinearDC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            upwind<Type>(mesh, faceFlux),
            gradScheme_
            (
                new fv::gaussGrad<Type>(mesh)
            )
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        skewCorrectedLinearDC
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    schemeData
                )
            )
        {}

        //- Construct from faceFlux and Istream
        skewCorrectedLinearDC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    schemeData
                )
            )
        {}


    // Member Functions

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        interpolate(const GeometricField<Type, fvPatchField, volMesh>&) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
