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

Description
    The scalar specialisations of the flux corrections, which have no meaning
    for a scalar "velocity" and fail as they do in every OpenFOAM ddt scheme.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

// The BDF ddt scheme is currently only implemented for OpenFOAM.com
#ifdef OPENFOAM_COM

#include "bdfDdtScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

bool bdfDdtSchemeWarned = false;


template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& U,
    const surfaceScalarField& Uf
)
{
    NotImplemented;
    return surfaceScalarField::null();
}


template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
)
{
    NotImplemented;
    return surfaceScalarField::null();
}


template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
)
{
    NotImplemented;
    return surfaceScalarField::null();
}


template<>
tmp<surfaceScalarField> bdfDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
)
{
    NotImplemented;
    return surfaceScalarField::null();
}

} // End namespace fv

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// makeFvDdtScheme is not used: it specialises the scalar flux corrections on
// the named class, which is done above
makeFvDdtTypeScheme(bdfDdtScheme, scalar)
makeFvDdtTypeScheme(bdfDdtScheme, vector)
makeFvDdtTypeScheme(bdfDdtScheme, sphericalTensor)
makeFvDdtTypeScheme(bdfDdtScheme, symmTensor)
makeFvDdtTypeScheme(bdfDdtScheme, tensor)

#endif

// ************************************************************************* //
