/*---------------------------------------------------------------------------*\
  XenoSim - Cardiovascular CFD & Multiphysics
  Component: BDFDdtScheme

  Author: David Kelly
  Copyright (C) 2026 David Kelly

Description
    The scalar specialisations of the flux corrections, which have no meaning
    for a scalar "velocity" and fail as they do in every OpenFOAM ddt scheme.

\*---------------------------------------------------------------------------*/

#include "BDFDdtScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

template<>
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& U,
    const surfaceScalarField& Uf
)
{
    NotImplemented;
    return surfaceScalarField::null();
}


template<>
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
)
{
    NotImplemented;
    return surfaceScalarField::null();
}


template<>
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtUfCorr
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
tmp<surfaceScalarField> BDFDdtScheme<scalar>::fvcDdtPhiCorr
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

// ************************************************************************* //
