/*---------------------------------------------------------------------------*\
  XenoSim - Cardiovascular CFD & Multiphysics
  Component: BDF4DdtScheme

  Author: David Kelly
  Copyright (C) 2026 David Kelly

Class
    Foam::fv::BDF4DdtScheme

Group
    grpFvDdtSchemes

Description
    Fourth-order backward differentiation formula ddt, using the current and
    four previous time-step values:

    \verbatim
        ddt(phi) = (25/12 phi^(n+1) - 4 phi^n + 3 phi^(n-1)
                  - 4/3 phi^(n-2) + 1/4 phi^(n-3))/dt
    \endverbatim

    Requires a constant time step and a static mesh. A fresh start ramps
    through Euler, BDF2 and BDF3 on its first three steps; see BDFDdtScheme
    for the start-up, and for how to supply the history so that it does not.

Usage
    \verbatim
    ddtSchemes
    {
        default         BDF4;
    }
    \endverbatim

    with \c libcirculationDdtSchemes loaded in system/controlDict.

SourceFiles
    BDF4DdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef BDF4DdtScheme_H
#define BDF4DdtScheme_H

#include "BDFDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class BDF4DdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BDF4DdtScheme
:
    public BDFDdtScheme<Type>
{
public:

    //- Runtime type information
    TypeName("BDF4");


    // Constructors

        //- Construct from mesh
        explicit BDF4DdtScheme(const fvMesh& mesh)
        :
            BDFDdtScheme<Type>(mesh, 4)
        {}

        //- Construct from mesh and Istream
        BDF4DdtScheme(const fvMesh& mesh, Istream& is)
        :
            BDFDdtScheme<Type>(mesh, is, 4)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
