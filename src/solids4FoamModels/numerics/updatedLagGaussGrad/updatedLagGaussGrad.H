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
    Foam::fv::updatedLagGaussGrad

Description
    Basic second-order gradient scheme using face-interpolation
    and Gauss' theorem.

SourceFiles
    updatedLagGaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef updatedLagGaussGrad_H
#define updatedLagGaussGrad_H

#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class updatedLagGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class updatedLagGaussGrad
:
    public fv::gradScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        updatedLagGaussGrad(const updatedLagGaussGrad&);

        //- Disallow default bitwise assignment
        void operator=(const updatedLagGaussGrad&);


public:

    //- Runtime type information
    TypeName("updatedLagGauss");


    // Constructors

        //- Construct from mesh
        updatedLagGaussGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linear<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        updatedLagGaussGrad(const fvMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(nullptr)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type> >
                    (
                        new linear<Type>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type> >
                    (
                        surfaceInterpolationScheme<Type>::New(mesh, is)
                    );
            }
        }


    // Member Functions

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        static
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > gradf
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const word& name
        );

        //- Return the gradient of the given field calculated
        //  using Gauss' theorem on the interpolated field
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;

        //- Return the BlockLduSystem corresponding to the implicit grad
        //  discretisation for coupledPressureDisplacementSolid
        tmp
        <
            BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
        > fvmGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        // Correct boundary conditions moved to base class
        // HJ, 14/Jun/2013
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarGaussGrad.H"

#ifdef NoRepository
#   include "updatedLagGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
