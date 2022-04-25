
#include "twoDPointCorrectorRBF.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    // Set marker to -1 if is on the "back" side
    void twoDPointCorrectorRBF::setMarker()
    {
        if ( !required() )
            return;

        // Change size of useablePointIDs_ and shadowPointIDs_
        useablePointIDs_.setSize( mesh_.nPoints() / 2 );
        shadowPointIDs_.setSize( mesh_.nPoints() / 2 );

        // Get reference to edges
        const edgeList & meshEdges = mesh_.edges();
        const pointField & points( mesh_.points() );

        const labelList & neIndices = normalEdgeIndices();
        const vector & pn = planeNormal();

        forAll( neIndices, edgeI )
        {
            const label & pStartInd = meshEdges[neIndices[edgeI]].start();
            const label & pEndInd = meshEdges[neIndices[edgeI]].end();
            const point & pStart = points[pStartInd];
            const point & pEnd = points[pEndInd];

            // calculate average point position
            const point A = 0.5 * (pStart + pEnd);

            // Calculate inner product with plane normal
            scalar pStartInner = ( pn & (pStart - A) );
            scalar pEndInner = ( pn & (pEnd - A) );

            if ( pStartInner > 0 && pEndInner < 0 )
            {
                pointMarker_[pEndInd] = -1;
                useablePointIDs_[edgeI] = pStartInd;
                shadowPointIDs_[edgeI] = pEndInd;
            }
            else
            if ( pEndInner > 0 && pStartInner < 0 )
            {
                pointMarker_[pStartInd] = -1;
                useablePointIDs_[edgeI] = pEndInd;
                shadowPointIDs_[edgeI] = pStartInd;
            }
            else
            {
                FatalErrorIn( "void twoDPointCorrectorRBF::setMarker()" )
                    << "Both points give back a negative value with the inner product. Programming error?"
                    << abort( FatalError );
            }
        }
    }

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    twoDPointCorrectorRBF::twoDPointCorrectorRBF( const polyMesh & mesh )
        :
        twoDPointCorrector( mesh ),
        mesh_( mesh ),
        pointMarker_( mesh.nPoints(), 0 ),
        useablePointIDs_( mesh.nPoints(), 0 ),
        shadowPointIDs_( mesh.nPoints(), 0 )
    {
        setMarker();
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    twoDPointCorrectorRBF::~twoDPointCorrectorRBF()
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    const labelList & twoDPointCorrectorRBF::marker() const
    {
        return pointMarker_;
    }

    void twoDPointCorrectorRBF::setShadowSide( vectorField & newpoints ) const
    {
        forAll( useablePointIDs_, ipoint )
        {
            newpoints[shadowPointIDs_[ipoint]] = newpoints[useablePointIDs_[ipoint]];
        }
    }
} // End namespace Foam
