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
    volumeSmoother

Description
    Volume smoothing, where each point is moved based on the neighbouring cell
    centres using weights based on the volumes of the cells.

SourceFiles
    volumeSmoother.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef volumeSmoother_H
#define volumeSmoother_H

#include "meshSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class volumeSmoother Declaration
\*---------------------------------------------------------------------------*/

class volumeSmoother
:
    public meshSmoother
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        volumeSmoother(const volumeSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const volumeSmoother&);

public:

    //- Runtime type information
    TypeName("volumeSmoother");

    // Static data members


    // Constructors

        //- Construct from dictionary
        volumeSmoother
        (
            const word& name,
            fvMesh& mesh
        );


    // Destructor

        virtual ~volumeSmoother();


    // Member Functions

        //- Smoothing function
        virtual scalar smooth(bool writeIters);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
