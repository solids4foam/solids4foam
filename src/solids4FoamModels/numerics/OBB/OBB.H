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
    Foam::OBB

Description
    Oriented bounding box represented by a centre, three orthonormal axes and
    the corresponding half-lengths.

    The box fitted to an arbitrary point field uses principal component
    analysis.  This generally gives a tighter box than an axis-aligned bound
    box, but is not guaranteed to be the minimum-volume oriented box.

    Box overlap is evaluated with the separating axis theorem.

References
    Eberly, D. H. (2008). 3D Game Engine Design (2nd ed.). Morgan Kaufmann.

    Ericson, C. (2005). Real-Time Collision Detection. Morgan Kaufmann.

    Gottschalk, S. (2000). Collision Queries Using Oriented Bounding Boxes.
    PhD thesis, University of North Carolina at Chapel Hill.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

SourceFiles
    OBB.C

\*---------------------------------------------------------------------------*/

#ifndef OBB_H
#define OBB_H

#include "boundBox.H"
#include "pointField.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OBB;

Ostream& operator<<(Ostream& os, const OBB& box);
Istream& operator>>(Istream& is, OBB& box);


/*---------------------------------------------------------------------------*\
                           Class OBB Declaration
\*---------------------------------------------------------------------------*/

class OBB
{
    // Private Data

        //- Centre of the box
        point centre_;

        //- Half-length along each local box axis
        vector halfLength_;

        //- Orthonormal local box axes stored as tensor rows
        tensor axes_;

        //- Angular tolerance used to identify parallel separating axes
        static const scalar parallelTolerance_;


    // Private Member Functions

        //- Fit a principal-component box to the supplied points
        void makeOBB(const pointField& points);

        //- Fit a minimum-area in-plane box and extrude it along faceNormal
        void makeOBB
        (
            const pointField& facePoints,
            const vector& faceNormal,
            const scalar normalPosExtent,
            const scalar normalNegExtent,
            const scalar scaleFactor
        );

        //- Calculate the point mean and symmetric covariance tensor
        static void calcCovariance
        (
            const pointField& points,
            point& mean,
            symmTensor& covariance
        );

        //- Expand local component bounds to include the supplied point
        static void addToLocalBounds
        (
            const vector& localPoint,
            vector& minLocal,
            vector& maxLocal
        );

        //- Orthonormalise the rows of axes_
        void orthonormaliseAxes();


public:

    // Static Data Members

        //- A box centred at the origin with VGREAT half-lengths
        static const OBB greatBox;


    // Constructors

        //- Construct a zero-sized box at the origin
        OBB();

        //- Construct an axis-aligned box from its centre and half-lengths
        OBB(const point& centre, const vector& halfLength);

        //- Construct from centre, half-lengths and local axes stored as rows
        OBB
        (
            const point& centre,
            const vector& halfLength,
            const tensor& axes
        );

        //- Construct from an axis-aligned bounding box
        OBB(const boundBox& box);

        //- Fit a principal-component box to the supplied points
        OBB(const pointField& points);

        //- Fit a principal-component box to a temporary point field
        OBB(const tmp<pointField>& points);

        //- Construct a face box with positive and negative normal extrusion
        OBB
        (
            const pointField& facePoints,
            const vector& faceNormal,
            const scalar normalPosExtent,
            const scalar normalNegExtent,
            const scalar scaleFactor = 1.0
        );

        //- Construct from Istream
        OBB(Istream& is);


    // Member Functions

        // Access

            //- Return the box centre
            const point& midpoint() const
            {
                return centre_;
            }

            //- Return the half-lengths along the local box axes
            const vector& ext() const
            {
                return halfLength_;
            }

            //- Return the local box axes stored as tensor rows
            const tensor& R() const
            {
                return axes_;
            }


        // Query

            //- Return true if this box overlaps the supplied box
            bool overlaps(const OBB& box) const;


    // Member Operators

        bool operator==(const OBB& box) const;

        bool operator!=(const OBB& box) const
        {
            return !operator==(box);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, OBB& box);
        friend Ostream& operator<<(Ostream& os, const OBB& box);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
