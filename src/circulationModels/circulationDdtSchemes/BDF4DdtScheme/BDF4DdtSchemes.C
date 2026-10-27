/*---------------------------------------------------------------------------*\
  XenoSim - Cardiovascular CFD & Multiphysics
  Component: BDF4DdtScheme

  Author: David Kelly
  Copyright (C) 2026 David Kelly
\*---------------------------------------------------------------------------*/

#include "BDF4DdtScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// makeFvDdtScheme is not used: it specialises the scalar flux corrections on
// the named class, and here they belong to the shared BDFDdtScheme base.
makeFvDdtTypeScheme(BDF4DdtScheme, scalar)
makeFvDdtTypeScheme(BDF4DdtScheme, vector)
makeFvDdtTypeScheme(BDF4DdtScheme, sphericalTensor)
makeFvDdtTypeScheme(BDF4DdtScheme, symmTensor)
makeFvDdtTypeScheme(BDF4DdtScheme, tensor)

// ************************************************************************* //
