/*---------------------------------------------------------------------------*\
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
    twoDPointCorrectorRBF

Description
    2D point corrector for RBF mesh motion.

Author
    David Blom.
    Philip Cardiff, UCD.

SourceFiles
    twoDPointCorrectorRBF.C

\*---------------------------------------------------------------------------*/

#ifndef twoDPointCorrectorRBF_H
#define twoDPointCorrectorRBF_H

#include "pointField.H"
#include "labelList.H"
#include "vector.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
*                       Class twoDPointCorrectorRBF Declaration
\*---------------------------------------------------------------------------*/

class twoDPointCorrectorRBF
:
    public twoDPointCorrector
{
    // Private data

        // - Reference to moving mesh
        const polyMesh & mesh_;

        // - Holder for marker value: -1 if on "back" side and 0 otherwise
        labelList pointMarker_;

        // - Point IDs marked with 0
        labelList useablePointIDs_;

        // - Point IDs marked with -1
        labelList shadowPointIDs_;


    // Private Member Functions

        // - Disallow default bitwise copy construct
        twoDPointCorrectorRBF( const twoDPointCorrectorRBF & );

        // - Disallow default bitwise assignment
        void operator=( const twoDPointCorrectorRBF & );

        void setMarker();

    public:

    // Constructors

        // - Construct from components
        twoDPointCorrectorRBF( const polyMesh & mesh );


    // Destructor

        virtual ~twoDPointCorrectorRBF();


    // Member Functions

        const labelList & marker() const;

        void setShadowSide(vectorField & newpoints) const;
};

} // End namespace Foam


#endif
