/*---------------------------------------------------------------------------*\
  XenoSim - Cardiovascular CFD & Multiphysics
  Component: BDF3DdtScheme

  Author: David Kelly
  Copyright (C) 2026 David Kelly
\*---------------------------------------------------------------------------*/

#include "BDF3DdtScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// makeFvDdtScheme is not used: it specialises the scalar flux corrections on
// the named class, and here they belong to the shared BDFDdtScheme base.
makeFvDdtTypeScheme(BDF3DdtScheme, scalar)
makeFvDdtTypeScheme(BDF3DdtScheme, vector)
makeFvDdtTypeScheme(BDF3DdtScheme, sphericalTensor)
makeFvDdtTypeScheme(BDF3DdtScheme, symmTensor)
makeFvDdtTypeScheme(BDF3DdtScheme, tensor)

// ************************************************************************* //
