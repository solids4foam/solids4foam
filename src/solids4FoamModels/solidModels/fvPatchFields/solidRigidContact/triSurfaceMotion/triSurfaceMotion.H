/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    triSurfaceMotion

Description
    Class in of motion for triSurfaces.

SourceFiles
    triSurfaceMotion.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMotion_H
#define triSurfaceMotion_H

#include "foamTime.H"
#include "typeInfo.H"
#include "standAlonePatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class triSurfaceMotion Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMotion
{
    // Private data

        //- Const reference to time
        const Time& runTime_;

        //- Copy of dictionary
        const dictionary dict_;

        //- Const reference to the standAlonePatch
        const standAlonePatch& patch_;

        //- Total displacement
        const vector displacement_;

        //- Displacement ramp time
        const scalar dispRampTime_;

        //- Maximum RPM
        const scalar rpm_;

        //- Rotation ramp time
        const scalar rotationRampTime_;

        //- Rotation axis
        const vector rotationAxis_;

        //- Current axis displacement
        vector currentAxisDisplacement_;

        //- Current angle
        scalar currentAngle_;

        //- Current RPM
        scalar currentRpm_;

        //- Original position of the patch points
        const pointField origPoints_;

        //- Original position of the face centres
        const vectorField origFaceCentres_;

        //- Old point displacements
        pointField oldPointDisp_;

        //- Old face displacements
        vectorField oldFaceDisp_;

        //- Initial rotation origin i.e. a point on the initial axis
        vector initialRotationOrigin_;

        // Current face displacement increments
        vectorField faceDisplacementIncrement_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        triSurfaceMotion(const triSurfaceMotion&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMotion&);


public:

    //- Runtime type information
    TypeName("triSurfaceMotion");


    // Constructors

        //- Construct from the time, dict and patch
        triSurfaceMotion
        (
            const Time& runTime,
            const dictionary& dict,
            const standAlonePatch& patch
        );

    // Destructor

        virtual ~triSurfaceMotion();


    // Member Functions

        //- Return the new patch points
        tmp<pointField> newPoints();

        //- Return the face displacement increments
        const vectorField& faceDisplacementIncrement() const
        {
            return faceDisplacementIncrement_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
