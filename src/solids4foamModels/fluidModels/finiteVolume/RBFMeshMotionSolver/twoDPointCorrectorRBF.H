
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

            void setShadowSide( vectorField & newpoints ) const;
    };
} // End namespace Foam


#endif
