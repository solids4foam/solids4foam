/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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
    newDynamicBodyFvMesh

Description
    Automatic motion of the mesh around a moving body. A direction,
    amplitude and frequency of translational motion and origin, axis,
    amplitude and frequency of rotational motion must be specified.

SourceFiles
    newDynamicBodyFvMesh.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef newDynamicBodyFvMesh_H
#define newDynamicBodyFvMesh_H

#include "dynamicFvMesh.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class newDynamicBodyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class newDynamicBodyFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary dynamicMeshCoeffs_;

        //- Moving mesh solver
        autoPtr<motionSolver> motionPtr_;

        //- Name of the patch which encloses the body
        word bodyPatchName_;

        //- Body patch ID
        label bodyPatchID_;

        //- Translational motion direction
        vector translationDirection_;

        //- Translational motion amplitude
        scalar translationAmplitude_;

        //- Translational motion frequency
        scalar translationFrequency_;

        //- Initial rotational motion origin
        vector initialRotationOrigin_;

        //- Rotational motion axis
        vector rotationAxis_;

        //- Rotational motion amplitude
        scalar rotationAmplitude_;

        //- Rotational motion frequency
        scalar rotationFrequency_;

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        newDynamicBodyFvMesh(const newDynamicBodyFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const newDynamicBodyFvMesh&);


public:

    //- Runtime type information
    TypeName("newDynamicBodyFvMesh");


    // Constructors

        //- Construct from IOobject
        newDynamicBodyFvMesh(const IOobject& io);


    // Destructor

        virtual ~newDynamicBodyFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
