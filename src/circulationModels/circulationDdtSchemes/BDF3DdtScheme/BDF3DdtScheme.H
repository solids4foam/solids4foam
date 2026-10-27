/*---------------------------------------------------------------------------*\
  XenoSim - Cardiovascular CFD & Multiphysics
  Component: BDF3DdtScheme

  Author: David Kelly
  Copyright (C) 2026 David Kelly

Class
    Foam::fv::BDF3DdtScheme

Group
    grpFvDdtSchemes

Description
    Third-order backward differentiation formula ddt, using the current and
    three previous time-step values:

    \verbatim
        ddt(phi) = (11/6 phi^(n+1) - 3 phi^n + 3/2 phi^(n-1) - 1/3 phi^(n-2))/dt
    \endverbatim

    Requires a constant time step and a static mesh. A fresh start ramps
    through Euler and BDF2 on its first two steps; see BDFDdtScheme for the
    start-up, and for how to supply the history so that it does not.

Usage
    \verbatim
    ddtSchemes
    {
        default         BDF3;
    }
    \endverbatim

    with \c libcirculationDdtSchemes loaded in system/controlDict.

SourceFiles
    BDF3DdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef BDF3DdtScheme_H
#define BDF3DdtScheme_H

#include "BDFDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class BDF3DdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BDF3DdtScheme
:
    public BDFDdtScheme<Type>
{
public:

    //- Runtime type information
    TypeName("BDF3");


    // Constructors

        //- Construct from mesh
        explicit BDF3DdtScheme(const fvMesh& mesh)
        :
            BDFDdtScheme<Type>(mesh, 3)
        {}

        //- Construct from mesh and Istream
        BDF3DdtScheme(const fvMesh& mesh, Istream& is)
        :
            BDFDdtScheme<Type>(mesh, is, 3)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
