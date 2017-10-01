/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    laplacianSmoother

Description
    Smoothing mesh using Laplacian smoothing. Needs a factor lambda to relax
    the smoothing in each iterations.

SourceFiles
    laplacianSmoother.C

Author
    Philip Cardiff UCD
    Zeljko Tukovic FSB/UCD

\*---------------------------------------------------------------------------*/

#ifndef laplacianSmoother_H
#define laplacianSmoother_H

#include "meshSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class laplacianSmoother Declaration
\*---------------------------------------------------------------------------*/

class laplacianSmoother
:
    public meshSmoother
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        laplacianSmoother(const laplacianSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const laplacianSmoother&);

public:

    //- Runtime type information
    TypeName("laplacianSmoother");

    // Static data members


    // Constructors

        //- Construct from dictionary
        laplacianSmoother
        (
            const word& name,
            fvMesh& mesh
        );


    // Destructor

        virtual ~laplacianSmoother();


    // Member Functions
        
        //- Smoothing function
        virtual scalar smooth(bool writeIters);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
