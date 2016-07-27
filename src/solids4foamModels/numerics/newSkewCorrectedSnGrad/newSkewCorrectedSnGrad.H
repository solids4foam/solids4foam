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
    Foam::fv::newSkewCorrectedSnGrad

Description
    Simple central-difference snGrad scheme with skewness and
    non-orthogonal correction.

SourceFiles
    newSkewCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef newSkewCorrectedSnGrad_H
#define newSkewCorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class newSkewCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class newSkewCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private data

        scalar limitCoeff_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const newSkewCorrectedSnGrad&);


public:

    //- Runtime type information
    TypeName("newSkewCorrected");


    // Constructors

        //- Construct from mesh
        newSkewCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        newSkewCorrectedSnGrad(const fvMesh& mesh, Istream& is)
        :
            snGradScheme<Type>(mesh),
            limitCoeff_(readScalar(is))
        {
            if (limitCoeff_ < 0 || limitCoeff_ > 1)
            {
                FatalIOErrorIn
                (
                    "newSkewCorrectedSnGrad(const fvMesh& mesh, Istream& is) : ",
                    is
                )   << "limitCoeff is specified as " << limitCoeff_
                    << " but should be >= 0 && <= 1"
                    << exit(FatalIOError);
            }
        }


    // Destructor

        virtual ~newSkewCorrectedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
//             return !this->mesh().orthogonal();
        }

        //- Return the explicit correction to the newSkewCorrectedSnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >
newSkewCorrectedSnGrad<scalar>::correction
(
    const GeometricField<scalar, fvPatchField, volMesh>&
) const;

template<>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> >
newSkewCorrectedSnGrad<vector>::correction
(
    const GeometricField<vector, fvPatchField, volMesh>&
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "newSkewCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
