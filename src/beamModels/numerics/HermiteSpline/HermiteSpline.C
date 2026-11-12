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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "HermiteSpline.H"
#include "scalarMatrices.H"
#include "demandDrivenData.H"
#include "LUscalarMatrix.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::HermiteSpline::calcParam()
{
    param_.setSize(points_.size());

    if (param_.size())
    {
        param_[0] = 0.0;

        scalarField segLen (segLengths());

        for (label i=1; i < param_.size(); i++)
        {
            param_[i] = param_[i-1] + segLen[i-1];
        }

        lineLength_ = param_.last();
    }
    else
    {
        lineLength_ = 0.0;
    }
}

void Foam::HermiteSpline::calcMidPoints() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (midPointsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::HermiteSpline::calcMidPoints() const"
        )
            << "Segment mid points already exist"
            << abort(FatalError);
    }

    midPointsPtr_ = new vectorField(nSegments(), vector::zero);
    vectorField& midPoints = *midPointsPtr_;

    forAll(midPoints, sI)
    {
        scalar segLen = param_[sI+1] - param_[sI];
        scalar midPointZeta = localParameter(sI, segLen/2);

        midPoints[sI] = position(sI, midPointZeta);
    }
}

void Foam::HermiteSpline::calcMidPointDerivatives() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (midPointDerivativesPtr_)
    {
        FatalErrorIn
        (
            "void Foam::HermiteSpline::calcMidPointDerivatives() const"
        )
            << "Segment mid point derivatives already exist"
            << abort(FatalError);
    }

    midPointDerivativesPtr_ = new vectorField(nSegments(), vector::zero);
    vectorField& midPointDerivatives = *midPointDerivativesPtr_;

    forAll(midPointDerivatives, sI)
    {
        scalar segLen = param_[sI+1] - param_[sI];
        scalar midPointZeta = localParameter(sI, segLen/2);

        midPointDerivatives[sI] = firstDerivative(sI, midPointZeta);
    }

    // midPointDerivatives /= mag(midPointDerivatives);
}

void Foam::HermiteSpline::calcMidPointCurvatures() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (midPointCurvaturesPtr_)
    {
        FatalErrorIn
        (
            "void Foam::HermiteSpline::calcMidPointCurvatures() const"
        )
            << "Segment mid point curvatures already exist"
            << abort(FatalError);
    }

    midPointCurvaturesPtr_ = new scalarField(nSegments(), 0);
    scalarField& midPointCurvatures = *midPointCurvaturesPtr_;

    forAll(midPointCurvatures, sI)
    {
        scalar segLen = param_[sI+1] - param_[sI];
        scalar midPointZeta = localParameter(sI, segLen/2);

        vector dRdt = paramFirstDerivative(sI, midPointZeta);
        vector d2Rdt2 = paramSecondDerivative(sI, midPointZeta);

        scalar tmp =
            (dRdt & dRdt)*(d2Rdt2 & d2Rdt2)
          - sqr(dRdt & d2Rdt2);

        if (tmp > SMALL)
        {
            midPointCurvatures[sI] =
                sqrt
                (
                    (
                        (dRdt & dRdt)*(d2Rdt2 & d2Rdt2)
                      - sqr(dRdt & d2Rdt2)
                    )
                   /(::pow((dRdt & dRdt), 3) + SMALL)
                );
        }
        else
        {
            midPointCurvatures[sI] = 0;
        }
    }
}

void Foam::HermiteSpline::calcDRdS() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (dRdSPtr_)
    {
        FatalErrorIn
        (
            "void Foam::HermiteSpline::calcDRdS() const"
        )
            << "Segment end point derivatives already exist"
            << abort(FatalError);
    }

    dRdSPtr_ = new vectorField(tangents_);
    // vectorField& dRdS = *dRdSPtr_;
    // dRdS /= mag(dRdS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HermiteSpline::HermiteSpline()
:
    points_(),
    tangents_(),
    c_(),
    midPointsPtr_(NULL),
    midPointDerivativesPtr_(NULL),
    midPointCurvaturesPtr_(NULL),
    dRdSPtr_(NULL),
    lineLength_(0.0),
    param_(0)
{}

Foam::HermiteSpline::HermiteSpline
(
    const pointField& ps,
    const pointField& ts
)
:
    points_(ps),
    tangents_(ts),
    c_(ps.size()-1, 1),
    midPointsPtr_(NULL),
    midPointDerivativesPtr_(NULL),
    midPointCurvaturesPtr_(NULL),
    dRdSPtr_(NULL),
    lineLength_(0.0),
    param_(0)
{
    tangents_ /= mag(tangents_);

    for (label i=0; i < c_.size(); i++)
    {
        c_[i] = mag(points_[i+1] - points_[i]);
    }

    calcParam();
}

Foam::HermiteSpline::HermiteSpline(Istream& s)
:
    points_(),
    tangents_(),
    c_(),
    midPointsPtr_(NULL),
    midPointDerivativesPtr_(NULL),
    midPointCurvaturesPtr_(NULL),
    dRdSPtr_(NULL),
    lineLength_(0.0),
    param_(0)
{
    points_ = vectorField(s);
    tangents_ = vectorField(s);

    tangents_ /= mag(tangents_);

    for (label i=0; i < c_.size(); i++)
    {
        c_[i] = mag(points_[i+1] - points_[i]);
    }

    calcParam();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HermiteSpline::~HermiteSpline()
{
    deleteDemandDrivenData(midPointsPtr_);
    deleteDemandDrivenData(midPointDerivativesPtr_);
    deleteDemandDrivenData(midPointCurvaturesPtr_);
    deleteDemandDrivenData(dRdSPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::HermiteSpline::points() const
{
    return points_;
}

const Foam::pointField& Foam::HermiteSpline::tangents() const
{
    return tangents_;
}

Foam::label Foam::HermiteSpline::nSegments() const
{
    return points_.size()-1;
}

Foam::label Foam::HermiteSpline::localParameter(scalar& lambda) const
{
    // check endpoints
    if (lambda < SMALL)
    {
        lambda = -1;
        return 0;
    }
    else if (lambda > lineLength_ - SMALL)
    {
        lambda = 1;
        return nSegments()-1;
    }

    // search table of cumulative distances to find which line-segment
    // we are on. Check the upper bound.

    label segmentI = 1;
    while (param_[segmentI] < lambda)
    {
        segmentI++;
    }
    segmentI--;   // we want the corresponding lower bound

    scalar arcLen = lambda - param_[segmentI];
    lambda = localParameter(segmentI, arcLen);

    return segmentI;
}

Foam::scalar Foam::HermiteSpline::localParameter
(
    const label segment,
    const scalar arcLen
) const
{
    scalar totalArcLength = param_[segment+1] - param_[segment];
    // scalar s = arcLength;

    if (arcLen > (totalArcLength-SMALL))
    {
        return 1;
    }
    else if (arcLen < SMALL)
    {
        return -1;
    }

    // Find zeta

    scalar residual = GREAT;
    scalar zeta0 = -1;
    scalar zeta1 =  1;
    scalar zeta2 =  0;
    do
    {
        scalar f0 = arcLength(segment, zeta0) - arcLen;
        scalar f1 = arcLength(segment, zeta1) - arcLen;

        zeta2 = zeta1 - f1*(zeta1 - zeta0)/(f1 - f0);

        residual = mag(zeta2 - zeta1)/2;

        zeta0 = zeta1;
        zeta1 = zeta2;
    }
    while(residual > 1e-4);

    return zeta2;
}

Foam::scalar Foam::HermiteSpline::arcLength
(
    const label segment,
    const scalar zeta
) const
{
    if (zeta < -1 + SMALL)
    {
        return 0;
    }

    // Dividing the arc length integral into further
    // subsegments
    label n(12*(zeta+1)/2);

    // n must be even number for Simpson's 1/3 rule to apply
    if ( (n % 2) > 0)
    {
        n += 1;
    }

    if (n < 2)
    {
        n = 2;
    }

    // Setting the length of step size between left end point
    // and zeta, i.e., [-1, zeta]
    const scalar h = (zeta + 1)/n;

    // Composite Simpson's 1/3 rule
    scalar arcLen = 0;
    for(label i = 1; i <= (n/2); i++)
    {
        arcLen +=
          h
         *(
              mag(paramFirstDerivative(segment, -1 + (2*i - 2)*h))
            + 4*mag(paramFirstDerivative(segment, -1 + (2*i - 1)*h))
            + mag(paramFirstDerivative(segment, -1 + (2*i)*h))
          )/3;
    }

    return arcLen;
}

Foam::point Foam::HermiteSpline::position(const scalar mu) const
{
    // check endpoints
    if (mu < SMALL)
    {
        return points_.first();
    }
    else if (mu > lineLength_ - SMALL)
    {
        return points_.last();
    }

    scalar zeta = mu;
    label segment = localParameter(zeta);
    return position(segment, zeta);
}


Foam::point Foam::HermiteSpline::position
(
    const label segment,
    const scalar zeta
) const
{
    // out-of-bounds
    if (segment < 0)
    {
        return points_.first();
    }
    else if (segment > nSegments())
    {
        return points_.last();
    }

    const point& p0 = points()[segment];
    const point& p1 = points()[segment+1];
    const point& t0 = tangents()[segment];
    const point& t1 = tangents()[segment+1];

    // special cases - no calculation needed
    if (zeta <= -1.0)
    {
        return p0;
    }
    else if (zeta >= 1.0)
    {
        return p1;
    }
    else
    {
        //Hermite polynomial interpolation
        scalar H0 = (2.0 + zeta)*sqr(1.0 - zeta)/4;
        scalar H1 = (2.0 - zeta)*sqr(1.0 + zeta)/4;
        scalar Ht0 =  (1.0 + zeta)*sqr(1.0 - zeta)/4;
        scalar Ht1 = -(1.0 - zeta)*sqr(1.0 + zeta)/4;

        vector curR = p0*H0 + p1*H1
          + c_[segment]*(Ht0*t0 + Ht1*t1)/2;

        return curR;
    }
}

Foam::point Foam::HermiteSpline::firstDerivative(const scalar mu) const
{
    // check endpoints
    if (mu < SMALL)
    {
        return tangents_.first();
    }
    else if (mu > lineLength_ - SMALL)
    {
        return tangents_.last();
    }

    scalar zeta = mu;
    label segment = localParameter(zeta);
    return firstDerivative(segment, zeta);
}

Foam::point Foam::HermiteSpline::firstDerivative
(
    const label segment,
    const scalar zeta
) const
{
    // out-of-bounds
    if (segment < 0)
    {
        return tangents_.first();
    }
    else if (segment > nSegments())
    {
        return tangents_.last();
    }

    const point& t0 = tangents()[segment];
    const point& t1 = tangents()[segment+1];

    // special cases - no calculation needed
    if (zeta <= -1)
    {
        return t0;
    }
    else if (zeta >= 1.0)
    {
        return t1;
    }
    else
    {
        scalar dH0dZeta = 3*(sqr(zeta) - 1)/4;
        scalar dH1dZeta = -3*(sqr(zeta) - 1)/4;
        scalar dHt0dZeta = (3*sqr(zeta) - 2*zeta - 1)/4;
        scalar dHt1dZeta = (3*sqr(zeta) + 2*zeta - 1)/4;

        const point& p0 = points()[segment];
        const point& p1 = points()[segment+1];

        vector dRdS = p0*dH0dZeta + p1*dH1dZeta
          + c_[segment]*(t0*dHt0dZeta + t1*dHt1dZeta)/2;

        dRdS /= mag(dRdS);

        return dRdS;
    }
}

Foam::point Foam::HermiteSpline::paramFirstDerivative
(
    const label segment,
    const scalar zeta
) const
{
    const point& t0 = tangents()[segment];
    const point& t1 = tangents()[segment+1];

    scalar dH0dZeta = 3*(sqr(zeta) - 1)/4;
    scalar dH1dZeta = -3*(sqr(zeta) - 1)/4;
    scalar dHt0dZeta = (3*sqr(zeta) - 2*zeta - 1)/4;
    scalar dHt1dZeta = (3*sqr(zeta) + 2*zeta - 1)/4;

    const point& p0 = points()[segment];
    const point& p1 = points()[segment+1];

    vector dRdZeta = p0*dH0dZeta + p1*dH1dZeta
      + c_[segment]*(t0*dHt0dZeta + t1*dHt1dZeta)/2;

    return dRdZeta;
}

Foam::point Foam::HermiteSpline::paramSecondDerivative
(
    const label segment,
    const scalar zeta
) const
{
    const point& t0 = tangents()[segment];
    const point& t1 = tangents()[segment+1];

    scalar d2H0dZeta2 =  3*zeta/2;
    scalar d2H1dZeta2 = -3*zeta/2;
    scalar d2Ht0dZeta2 = (3*zeta - 1)/2;
    scalar d2Ht1dZeta2 = (3*zeta + 1)/2;

    const point& p0 = points()[segment];
    const point& p1 = points()[segment+1];

    vector d2RdZeta2 = p0*d2H0dZeta2 + p1*d2H1dZeta2
      + c_[segment]*(t0*d2Ht0dZeta2 + t1*d2Ht1dZeta2)/2;

    return d2RdZeta2;
}

// Foam::scalar Foam::HermiteSpline::Jacobian
// (
//     const label segment,
//     const scalar zeta
// ) const
// {
//     scalar dH0dZeta = 3*(sqr(zeta) - 1)/4;
//     scalar dH1dZeta = -3*(sqr(zeta) - 1)/4;
//     scalar dHt0dZeta = (3*sqr(zeta) - 2*zeta - 1)/4;
//     scalar dHt1dZeta = (3*sqr(zeta) + 2*zeta - 1)/4;

//     const point& p0 = points()[segment];
//     const point& p1 = points()[segment+1];
//     const point& t0 = tangents()[segment];
//     const point& t1 = tangents()[segment+1];

//     vector dRdZeta = p0*dH0dZeta + p1*dH1dZeta
//       + c_[segment]*(t0*dHt0dZeta + t1*dHt1dZeta)/2;

//     return mag(dRdZeta);
// }

Foam::scalar Foam::HermiteSpline::length() const
{
    return lineLength_;
}

Foam::scalar Foam::HermiteSpline::length(const label segI, const scalar zeta) const
{
    return (param_[segI] + arcLength(segI, zeta));
}

const Foam::vectorField& Foam::HermiteSpline::midPoints() const
{
    if (!midPointsPtr_)
    {
        calcMidPoints();
    }

    return *midPointsPtr_;
}

const Foam::vectorField& Foam::HermiteSpline::dRdS() const
{
    if (!dRdSPtr_)
    {
        calcDRdS();
    }

    return *dRdSPtr_;
}

const Foam::vectorField& Foam::HermiteSpline::midPointDerivatives() const
{
    if (!midPointDerivativesPtr_)
    {
        calcMidPointDerivatives();
    }

    return *midPointDerivativesPtr_;
}

const Foam::scalarField& Foam::HermiteSpline::midPointCurvatures() const
{
    if (!midPointCurvaturesPtr_)
    {
        calcMidPointCurvatures();
    }

    return *midPointCurvaturesPtr_;
}

Foam::tmp<Foam::vectorField> Foam::HermiteSpline::firstDerivativeParam
(
    const labelScalarList& points
) const
{
    tmp<vectorField> tFirstDerivative
    (
        new vectorField(points.size(), vector::zero)
    );

    forAll(tFirstDerivative(), pI)
    {
        tFirstDerivative.ref()[pI] =
            paramFirstDerivative
            (
                points[pI].first(),
                points[pI].second()
            );
    }

    return tFirstDerivative;
}

Foam::tmp<Foam::vectorField> Foam::HermiteSpline::secondDerivativeParam
(
    const labelScalarList& points
) const
{
    tmp<vectorField> tSecondDerivative
    (
        new vectorField(points.size(), vector::zero)
    );

    forAll(tSecondDerivative(), pI)
    {
        tSecondDerivative.ref()[pI] =
            paramSecondDerivative
            (
                points[pI].first(),
                points[pI].second()
            );
    }

    return tSecondDerivative;
}


Foam::labelScalar Foam::HermiteSpline::nearestPoint
(
    const label segI,
    const vector& p
) const
{
    vector np(0, 0, 0);

    scalar mu0 = -1;
    scalar f0 = (tangents_[segI] & (p - points_[segI]));

    scalar mu1 = 1;
    scalar f1 = (tangents_[segI+1] & (p - points_[segI+1]));

    scalar mu = -2;

    if ((f0*f1) > SMALL)
    {
        FatalErrorIn
	    (
	        "Foam::HermiteSpline::nearestPoint(...) const"
	    )
	    << "Nearest point is out of segment"
	    << abort(FatalError);
    }
    else if (mag(f0) < SMALL)
    {
        np = points_[segI];
	    mu = mu0;
    }
    else if (mag(f1) < SMALL)
    {
        np = points_[segI + 1];
	    mu = mu1;
    }
    else // find root using Newton-Raphson
    {
        scalar err = GREAT;
        label nIter = 0;

        do
        {
            scalar df1 = (f1 - f0)/(mu1 - mu0);
            mu = mu1 - (f1/df1);
            scalar f =
            (
                firstDerivative(segI, mu)
              & (p - position(segI, mu))
            );

            if (f0*f > 0)
            {
                err = mag(mu - mu0);
                mu0 = mu;
                f0 = f;
            }
            else
            {
                err = mag(mu - mu1);
                mu1 = mu;
                f1 = f;
            }
        }

        while( (err > 1e-3) && (++nIter < 100) );

        // When nIter = 100, the while loop is exhausted, but still err > 1e-3,
        // should I need to add a Fatal Error here, or simply accept the last refined
        // value of mu? nearestPoint is used for line contact search - need to revisit this.

	    // if
	    // (
	    //     (mu < -1)
	    //     && (mu > 1)
	    // )
	    // {
        //     Info << "Nearest point is not found" << endl;
	    // }

        // SB changed this - 01 May 2024 (mu can never be both < -1 and > 1 at the same time)
        if
	    (
	        (mu < -1)
	        || (mu > 1)
	    )
	    {
            Info << "Nearest point is not found" << endl;
	    }

        np = position(segI, mu);
    }

    return labelScalar(segI, mu);
}

Foam::labelScalar Foam::HermiteSpline::findNearestPoint
(
    const label segI,
    const vector& p,
    const label nNeighbours
) const
{
    if (nNeighbours < 0)
    {
        FatalErrorIn
	    (
	     "Foam::HermiteSpline::findNearestPoint(...) const"
	    )
	    << "Number of neighbour segments must be >= 0"
	    << abort(FatalError);
    }

    labelScalar np(-1, -2);

    {
        // First segment
        label fSegI = segI - nNeighbours;
	    if (fSegI < 0)
	    {
	        fSegI = 0;
	    }

        // Last segment
        label lSegI = segI + nNeighbours;
        if (lSegI > (nSegments() - 1))
        {
                lSegI = nSegments() - 1;
        }

        for (label sI = fSegI; sI <= lSegI; sI++)
        {
            scalar testValue =
            (
                dRdS()[sI]
                & (p - points()[sI])
            )
            *(
                dRdS()[sI + 1]
                & (p - points()[sI + 1])
            );

            if (testValue <= 0)
            {
                np = nearestPoint(sI, p);
                break;
            }
        }

        if (np.first() != -1)
        {
            return np;
        }
        else
        {
            FatalErrorIn
            (
                "Foam::HermiteSpline::findNearestPoint(...) const"
            )
            << "Problem with finding nearest point"
            << abort(FatalError);
        }
    }

    return np;
}

Foam::Tuple2<Foam::scalar, Foam::scalar>
Foam::HermiteSpline::checkPointContact
(
    const label segI,
    const HermiteSpline& neiSpline,
    const label neiSegI,
    const scalar lowerContactAngleLimit
) const
{
    scalar zetaC = -2;
    scalar neiZetaC = -2;

    // Calculating potential contact angle
    vector tan0 = firstDerivative(segI, 0);
    tan0 /= mag(tan0) + SMALL;

    vector tan1 = neiSpline.firstDerivative(neiSegI, 0);
    tan1 /= mag(tan1) + SMALL;

    scalar contactAngle = mag(::acos(tan0 & tan1)*180.0/M_PI);

    if (contactAngle > lowerContactAngleLimit)
    {
        // Find initial solution point contact pair [zetaC,neiZetaC] using the
        // assumption of linear interpolation between points p0 and p1.
        // See Eq. 10 of the paper, Wriggers, P. and Zavarise, G., 1997. On contact between
        // three‐dimensional beams undergoing large deflections.
        // Communications in numerical methods in engineering, 13(6), pp.429-438.

        const vector& p0 = points()[segI];
        const vector& p1 = points()[segI+1];

        const vector& barP0 = neiSpline.points()[neiSegI];
        const vector& barP1 = neiSpline.points()[neiSegI+1];

        vector b = p1 + p0;
        vector barB = barP1 + barP0;

        vector t = p1 - p0;
        vector barT = barP1 - barP0;

        zetaC =
           -(
                (barB - b)
              & (
                    (barT*(barT & t) - t*(barT & barT))
                   /((barT & barT)*(t & t) - sqr(barT & t))
                )
            );

        neiZetaC =
            (
                (barB - b)
              & (
                    (t*(barT & t) - barT*(t & t))
                   /((barT & barT)*(t & t) - sqr(barT & t))
                )
            );

        // Then refine zetaC and neiZetaC values for Hermite spline interpolation
        // using Newton-Raphson technique.
        // Note - Revisit the if loop - possible reason why numerical oscillations
        // occur when point contact exactly occurs at a face
        if
        (
            (zetaC >= -1.05) && (zetaC <= 1.05)
         && (neiZetaC >= -1.05) && (neiZetaC <= 1.05)
        )
        {
            scalar residual = GREAT;
            label nIter = 0;
            do
            {
                nIter++;

                vector Rc = position(segI, zetaC);
                vector dRc = paramFirstDerivative(segI, zetaC);
                vector d2Rc = paramSecondDerivative(segI, zetaC);

                vector neiRc = neiSpline.position(neiSegI, neiZetaC);
                vector neiDRc =
                    neiSpline.paramFirstDerivative(neiSegI, neiZetaC);
                vector neiD2Rc =
                    neiSpline.paramSecondDerivative(neiSegI, neiZetaC);

                // Check Eq. 6 of the reference paper
                scalarSquareMatrix M(2, 0.0);
                M[0][0] = -(dRc & dRc) + ((neiRc - Rc) & d2Rc);
                M[0][1] =  (dRc & neiDRc);
                M[1][0] = -(dRc & neiDRc);
                M[1][1] =  (neiDRc & neiDRc) + ((neiRc - Rc) & neiD2Rc);

                scalarField rhs(2, 0);
                rhs[0] = -((neiRc - Rc) & dRc);
                rhs[1] = -((neiRc - Rc) & neiDRc);

                //scalarSquareMatrix invM = M.LUinvert();
		        LUscalarMatrix LUOfM(M);
                scalarSquareMatrix invM;
		        LUOfM.inv(invM);

                scalarField DZeta(2, 0);
                for (label i = 0; i < 2; i++)
                {
                    for (label j = 0; j< 2; j++)
                    {
                        DZeta[i] += invM[i][j]*rhs[j];
                    }
                }

                zetaC += DZeta[0];
                neiZetaC += DZeta[1];

                residual = max(mag(DZeta[0]/2), mag(DZeta[1]/2));

            }
            while(residual > 1e-4 && nIter < 100);
        }

    }
    // Comment SB - need to revisit 01 May 2024:
    // Can an else loop be added here to make it more clear
    // that if contactangle < lowerContactAngle limit then
    // multiple point contact points might occur, and single
    // point contact can never be found.

    // Also nIter=100 in CheckPointContact loop might not
    // be enough for angles (between 10-15 deg especially)
    // around the lowerContactAngle limit.
    // Beaware of this, ask ZT.

    return Tuple2<scalar, scalar>(zetaC, neiZetaC);
}

Foam::scalar Foam::HermiteSpline::segLength(const label segI) const
{
    return arcLength(segI, 1);
}

Foam::tmp<Foam::scalarField> Foam::HermiteSpline::segLengths() const
{
    tmp<scalarField> tSegLengths
    (
        new scalarField(nSegments(), 0)
    );

    forAll(tSegLengths(), segI)
    {
        tSegLengths.ref()[segI] = arcLength(segI, 1);
    }

    return tSegLengths;
}

Foam::tmp<Foam::scalarField> Foam::HermiteSpline::midPointParameters() const
{
    tmp<scalarField> tMidParameters
    (
        new scalarField(nSegments(), 0)
    );

    const scalarField segLen (segLengths());

    forAll(tMidParameters(), segI)
    {
        tMidParameters.ref()[segI] = param_[segI] + segLen[segI]/2;
    }

    return tMidParameters;
}

Foam::tmp<Foam::vectorField> Foam::HermiteSpline::segmentToPointInterpolate
(
    const vectorField& midPhi
)
const
{
    tmp<vectorField> tResult
    (
        new vectorField(points_.size(), vector::zero)
    );
     vectorField& result = tResult.ref();

    // Calculate tangent
    vectorField Dt(result.size(), vector::zero);
    Dt[0] = midPhi[1] - midPhi[0];
    Dt[0] /= mag(Dt[0]) + SMALL;
    Dt[result.size() - 1] =
        midPhi[result.size() - 1] - midPhi[result.size() - 2];
    Dt[result.size()-1] /= mag(Dt[result.size()-1]) + SMALL;
    for (label i=1; i<(Dt.size() - 1); i++)
    {
        Dt[i] = midPhi[i] - midPhi[i - 1];
        Dt[i] /= mag(Dt[i]) + SMALL;
    }

    // First and last points are fixed
    result[0] = midPhi[0];
    result[result.size() - 1] = midPhi[midPhi.size() - 1];
    for (label i = 1; i < (midPhi.size() - 1); i++)
    {
        scalar segLen = param_[i] - param_[i-1];
        scalar zeta = localParameter(i - 1, segLen/2);

        scalar H0 = (2.0 + zeta)*sqr(1.0 - zeta)/4;
        scalar H1 = (2.0 - zeta)*sqr(1.0 + zeta)/4;
        scalar Ht0 =  (1.0 + zeta)*sqr(1.0 - zeta)/4;
        scalar Ht1 = -(1.0 - zeta)*sqr(1.0 + zeta)/4;

        result[i] = midPhi[i] - result[i - 1]*H0
          - c_[i - 1]*(Dt[i-1]*Ht0 + Dt[i]*Ht1)/2;
        result[i] /= H1;
    }

    return tResult;
}

void Foam::HermiteSpline::movePoints
(
    const pointField& newPoints,
    const pointField& newTangents
)
{
    if (newPoints.size() != points_.size())
    {
        FatalErrorIn
        (
            "Foam::HermiteSpline::movePoints(...) const"
        )
            << "Number of new points is not equal to number of old points "
            << newPoints.size() << " != " << points_.size()
            << abort(FatalError);
    }

    points_ = newPoints;
    tangents_ = newTangents;
    tangents_ /= mag(tangents_);

    for (label i=0; i < c_.size(); i++)
    {
        c_[i] = mag(points_[i+1] - points_[i]);
    }

    calcParam();

    deleteDemandDrivenData(midPointsPtr_);
    deleteDemandDrivenData(midPointDerivativesPtr_);
    deleteDemandDrivenData(dRdSPtr_);
}

Foam::point Foam::HermiteSpline::position
(
    const vector& p0,
    const vector& t0,
    const vector& p1,
    const vector& t1,
    const scalar zeta
)
{
    // special cases - no calculation needed
    if (zeta <= -1.0)
    {
        return p0;
    }
    else if (zeta >= 1.0)
    {
        return p1;
    }
    else
    {
        //Hermite polynomial interpolation
        scalar H0 = (2.0 + zeta)*sqr(1.0 - zeta)/4;
        scalar H1 = (2.0 - zeta)*sqr(1.0 + zeta)/4;
        scalar Ht0 =  (1.0 + zeta)*sqr(1.0 - zeta)/4;
        scalar Ht1 = -(1.0 - zeta)*sqr(1.0 + zeta)/4;

        scalar c = mag(p1 - p0);

        vector curR = p0*H0 + p1*H1
          + c*(Ht0*t0 + Ht1*t1)/2;

        return curR;
    }
}


// ************************************************************************* //
