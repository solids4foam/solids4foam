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
    Foam::HermiteSpline

Description
    Hermite splines are cubic polynomials to interpolate between points
    using the positions of the points plus their tangents (derivatives).
    The basic form involves defining each segment of the spline by the start
    and end points (positions) and the tangents at these points. Here,
    for a line segment, the local parameter zeta [-1, 1] is used. zeta = -1 for left
    point and zeta = 1 for right end point.

    Advantages:
    1. Control: over the slope of curve at control points, advantegeous
        for ensuring smooth transitions (C-1 derivatives).
    2. Flexibility: Useful in applications where the derivative of the curve is known or
        needs to be controlled, such as in motion control paths in animation and robotics.
    3. Local Control: Adjusting the position or tangent at one point does not affect the entire
        curve, just the segment connected to that point.

Author
    Zeljko Tukovic, FSB Zagreb
    Seevani Bali, UCD.

SourceFiles
    HermiteSpline.C

\*---------------------------------------------------------------------------*/

#ifndef HermiteSpline_H
#define HermiteSpline_H

#include "pointField.H"
#include "scalarList.H"
#include "labelScalar.H"
#include "LUscalarMatrix.H"
#include "scalarMatrices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HermiteSpline Declaration
\*---------------------------------------------------------------------------*/

class HermiteSpline
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        HermiteSpline(const HermiteSpline&);

        //- Disallow default bitwise assignment
        void operator=(const HermiteSpline&);

protected:

    // Protected data

        //- The end points of each segments
        pointField points_;

        //- The end points tangents
        pointField tangents_;

        //- Hermite parameter
        scalarField c_;

        //- The mid points of segments
        mutable vectorField* midPointsPtr_;

        //- The segment mid point derivatives
        mutable vectorField* midPointDerivativesPtr_;

        //- The segment mid point curvatures
        mutable scalarField* midPointCurvaturesPtr_;

        //- The derivative at end points
        mutable vectorField* dRdSPtr_;

        //- The real line length
        scalar lineLength_;

        //- The rational (0-1) cumulative parameter value for each point
        scalarList param_;

    // Protected Member Functions

        //- Precalculate the rational cumulative parameter value
        //  and the line-length
        void calcParam();

        //- Return the line segment and the local parameter [0..1]
        //  corresponding to the global lambda [0..1]
        label localParameter(scalar& lambda) const;

        scalar localParameter
        (
            const label segment,
            const scalar arcLength
        ) const;

        //- Calculate the length of a line segment between interval [-1, zeta]
        // by integrating the first order derivatives of position vector using
        // composite Simpson's 1/3 rule
        // (https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_1/3_rule)
        scalar arcLength
        (
            const label segment,
            const scalar zeta
        ) const;

        //- Calculate segment mid points (zeta = 0)
        void calcMidPoints() const;

        //- Calculate segment mid point derivatives
        void calcMidPointDerivatives() const;

        //- Calculate segment mid point curvatures
        void calcMidPointCurvatures() const;

        //- Calculate segment end point derivatives
        void calcDRdS() const;

public:

    // Constructors

        //-
        HermiteSpline();

        //- Construct from components
        HermiteSpline
        (
            const pointField& points,
            const pointField& tangents
        );

        //- Construct from components
        HermiteSpline
        (
            Istream& stream
        );

    // Destructor

        virtual ~HermiteSpline();

    // Member Functions

        //- Return const-access to the control-points
        const pointField& points() const;

        //- Return const-access to the tangents
        const pointField& tangents() const;

        //- Return the number of line segments
        label nSegments() const;

        //- Return the point position corresponding to the curve parameter
        point position(const scalar) const;

        //- Return the point position corresponding to the local parameter
        point position(const label segment, const scalar) const;

        //- Return the first derivative corresponding to the curve parameter
        //  (cord length)
        point firstDerivative(const scalar) const;

        //- Return the first derivative corresponding to the local parameter
        point firstDerivative(const label segment, const scalar) const;

        //- Return the first derivative corresponding to the local parameter
        point paramFirstDerivative(const label segment, const scalar) const;

        //- Return the second derivative corresponding to the local parameter
        point paramSecondDerivative(const label segment, const scalar) const;

        // //- Return the Jacobian corresponding to the local parameter
        // scalar Jacobian(const label segment, const scalar) const;

        //- Return the length of the curve
        scalar length() const;

        //- Return the length of the curve
        scalar length(const label segI, const scalar zeta) const;

        //- Return midpoints of spline segments
        const vectorField& midPoints() const;

        //- Return segment midpoint derivatives
        const vectorField& midPointDerivatives() const;

        //- Return segment midpoint derivatives
        const scalarField& midPointCurvatures() const;

        //- Return first derivative for specified points
        //  (segment ID and parameter)
        tmp<vectorField> firstDerivativeParam
        (
            const labelScalarList& points
        ) const;

        //- Return second derivative for specified points
        //  (segment ID and parameter)
        tmp<vectorField> secondDerivativeParam
        (
            const labelScalarList& points
        ) const;

        //- Return tangents of spline segments
        const vectorField& dRdS() const;

        //- Nearest segment point for given point
        labelScalar nearestPoint
	    (
	        const label segI,
	        const vector& point
	    ) const;

        //- Nearest point for given point
        labelScalar findNearestPoint
	    (
	        const label segI,
	        const vector& point,
	        const label nNeighbours = 0
	    ) const;

        //- Check for point contact between two splines
        Tuple2<scalar, scalar> checkPointContact
        (
            const label segI,
            const HermiteSpline& neiSpline,
            const label neiSegI,
            const scalar contactAngleLimit
        ) const;

        //- Return segment lengths
        tmp<scalarField> segLengths() const;

        //- Return segment length
        scalar segLength(const label segI) const;

        //- Return midpoint paremeters (arc-length) of spline segments
        tmp<scalarField> midPointParameters() const;

        //- Interpolate field from segments to points
        tmp<vectorField> segmentToPointInterpolate
        (
            const vectorField& midPhi
        ) const;

        //- Move points
        void movePoints
        (
            const pointField& newPoints,
            const pointField& newTangents
        );

        //- Return position along segment for provided
        //  end points and tangents
        static point position
        (
            const vector& p0,
            const vector& t0,
            const vector& p1,
            const vector& t1,
            const scalar zeta
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
