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
    taubinSmoother

Description
    Smoothing mesh using Laplacian (Taubin) smoothing, with a shrink step
    followed by an inflation step.

    All points are smoothed, where boundary patch points are only allowed to
    slide along the boundary, and edge points are only allowed to slide along
    the edges.

    By default, the patches and edges are smoothed with the internal points,
    but optional switches allow them to be smoothed independently.

    Also, an edge point is not smoothed if the angle along the edge at that
    point is greater than the specified critical angle; by default it is 45
    degrees, but this can be optionally specified.

SourceFiles
    taubinSmoother.C

Author
    Philip Cardiff UCD
    Zeljko Tukovic FSB/UCD

\*---------------------------------------------------------------------------*/

#ifndef taubinSmoother_H
#define taubinSmoother_H

#include "meshSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class taubinSmoother Declaration
\*---------------------------------------------------------------------------*/

class taubinSmoother
:
    public meshSmoother
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        taubinSmoother(const taubinSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const taubinSmoother&);

public:

    //- Runtime type information
    TypeName("taubin");

    // Static data members


    // Constructors

        //- Construct from dictionary
        taubinSmoother
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~taubinSmoother();


    // Member Functions

        //- Evaluate different weighing options
        scalar weight(const word& weighingFunctionName, const vector& d);

        //- Smoothing function
        virtual scalar smooth();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
