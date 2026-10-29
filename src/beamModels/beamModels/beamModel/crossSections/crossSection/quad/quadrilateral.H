/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::quadrilateral

Description
    A quadrilateral primitive used to calculate face ...

SourceFiles
    quadrilateral.C

\*---------------------------------------------------------------------------*/

#ifndef quadrilateral_H
#define quadrilateral_H

#include "vector.H"
#include "tensor2D.H"
#include "Field.H"
#include "FixedList.H"
#include "pointField.H"
#include "face.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
// class Istream;
// class Ostream;

// Forward declaration of friend functions and operators

// class quadrilateral;

/*---------------------------------------------------------------------------*\
                           class quadrilateral Declaration
\*---------------------------------------------------------------------------*/

class quadrilateral
:
    // public face
    public FixedList<label, 4>
{
    // Private data

    // Private member functions

public:

    // Constructors

        //- Construct from four points
        quadrilateral
        (
            const label a,
            const label b,
            const label c,
            const label d
        );

        //- Construct from four points
        quadrilateral();

        //- Construct from four points
        quadrilateral(const labelList& pts);

        // //- Construct from Istream
        // inline quadrilateral(Istream&);

    // Member Functions

        //- Calculate average value at centroid of face
        template<class Type>
        Type integral(const Field<vector>&, const Field<Type>&) const;

        // Integrate nxm matrix field
        scalarRectangularMatrix integral
        (
            const Field<vector>& pts,
            const List<scalarRectangularMatrix>& A,
            const List<scalarRectangularMatrix>& CT
        ) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type integral
        (
            const Field<vector>&,
            const Field<Type>&,
            const Field<scalar>&
        ) const;

        // Access

            // //- Return first vertex
            // inline const vector& a() const;

            // //- Return second vertex
            // inline const vector& b() const;

            // //- Return third vertex
            // inline const vector& c() const;

            // //- Return fourth vertex
            // inline const vector& d() const;

        // Properties

            //- Return the points corresponding to this face
            inline pointField points(const pointField& points) const;

            // //- Return centre (centroid)
            // inline Point centre() const;

            //- Return scalar magnitude
            scalar area(const pointField& pts) const;

            // //- Return vector normal
            // inline vector normal() const;

            //- Evaluate shape function at a given point
            scalar Ni(label i, scalar xi, scalar eta) const;

            //- Evaluate shape function at a given point
            scalar dNidEta(label i, scalar xi, scalar eta) const;

            //- Evaluate shape function at a given point
            scalar dNidXi(label i, scalar xi, scalar eta) const;

            //- Evaluate shape function at a given point
            scalar detJacob
            (
                scalar xi,
                scalar eta,
                const pointField& points
            ) const;

    // IOstream operators

        // friend Istream& operator>>
        // (
        //     Istream&,
        //     quadrilateral&
        // );

        // friend Ostream& operator<<
        // (
        //     Ostream&,
        //     const quadrilateral&
        // );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "quadTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
