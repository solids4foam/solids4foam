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

\*---------------------------------------------------------------------------*/

#include "OBB.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::OBB::parallelTolerance_ = 1e-6;

const Foam::OBB Foam::OBB::greatBox
(
    Foam::point::zero,
    Foam::vector(Foam::VGREAT, Foam::VGREAT, Foam::VGREAT)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::OBB::calcCovariance
(
    const pointField& points,
    point& mean,
    symmTensor& covariance
)
{
    if (points.empty())
    {
        FatalErrorInFunction
            << "Cannot fit an OBB to an empty point field"
            << abort(FatalError);
    }

    mean = point::zero;

    forAll(points, pointI)
    {
        mean += points[pointI];
    }

    mean /= points.size();

    covariance = symmTensor::zero;

    forAll(points, pointI)
    {
        const vector d(points[pointI] - mean);

        covariance.xx() += d.x()*d.x();
        covariance.xy() += d.x()*d.y();
        covariance.xz() += d.x()*d.z();
        covariance.yy() += d.y()*d.y();
        covariance.yz() += d.y()*d.z();
        covariance.zz() += d.z()*d.z();
    }

    // Scaling does not change the eigenvectors. Avoid division for one point.
    if (points.size() > 1)
    {
        covariance /= scalar(points.size() - 1);
    }
}


void Foam::OBB::addToLocalBounds
(
    const vector& localPoint,
    vector& minLocal,
    vector& maxLocal
)
{
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        minLocal[cmpt] = Foam::min(minLocal[cmpt], localPoint[cmpt]);
        maxLocal[cmpt] = Foam::max(maxLocal[cmpt], localPoint[cmpt]);
    }
}


void Foam::OBB::orthonormaliseAxes()
{
    vector axis0(axes_.x());
    vector axis1(axes_.y());

    if (mag(axis0) < SMALL)
    {
        axis0 = vector(1, 0, 0);
    }
    else
    {
        axis0 /= mag(axis0);
    }

    axis1 -= (axis0 & axis1)*axis0;

    if (mag(axis1) < SMALL)
    {
        // Choose the Cartesian direction least aligned with axis0.
        if (mag(axis0.x()) <= mag(axis0.y()) && mag(axis0.x()) <= mag(axis0.z()))
        {
            axis1 = vector(1, 0, 0);
        }
        else if (mag(axis0.y()) <= mag(axis0.z()))
        {
            axis1 = vector(0, 1, 0);
        }
        else
        {
            axis1 = vector(0, 0, 1);
        }

        axis1 -= (axis0 & axis1)*axis0;
    }

    axis1 /= mag(axis1);

    vector axis2(axis0 ^ axis1);
    axis2 /= mag(axis2);

    axes_ = tensor(axis0, axis1, axis2);
}


void Foam::OBB::makeOBB(const pointField& points)
{
    point mean(point::zero);
    symmTensor covariance(symmTensor::zero);
    calcCovariance(points, mean, covariance);

    axes_ = eigenVectors(covariance);
    orthonormaliseAxes();

    vector minLocal(VGREAT, VGREAT, VGREAT);
    vector maxLocal(-VGREAT, -VGREAT, -VGREAT);

    forAll(points, pointI)
    {
        addToLocalBounds
        (
            axes_ & (points[pointI] - mean),
            minLocal,
            maxLocal
        );
    }

    const vector localCentre(0.5*(minLocal + maxLocal));
    centre_ = mean + (localCentre & axes_);
    halfLength_ = 0.5*(maxLocal - minLocal);
}


void Foam::OBB::makeOBB
(
    const pointField& facePoints,
    const vector& faceNormal,
    const scalar normalPosExtent,
    const scalar normalNegExtent,
    const scalar scaleFactor
)
{
    if (facePoints.size() < 3)
    {
        FatalErrorInFunction
            << "At least three face points are required; supplied "
            << facePoints.size() << abort(FatalError);
    }

    if (mag(faceNormal) < SMALL)
    {
        FatalErrorInFunction
            << "The face normal has zero magnitude" << abort(FatalError);
    }

    if (normalPosExtent < 0 || normalNegExtent < 0 || scaleFactor <= 0)
    {
        FatalErrorInFunction
            << "The normal extents must be non-negative and scaleFactor must "
            << "be positive" << abort(FatalError);
    }

    const vector normal(faceNormal/mag(faceNormal));
    scalar minArea = VGREAT;
    vector bestAxis1(vector::zero);
    vector bestAxis2(vector::zero);
    point bestCentre(point::zero);

    forAll(facePoints, pointI)
    {
        const label previousI = pointI ? pointI - 1 : facePoints.size() - 1;
        vector axis1(facePoints[pointI] - facePoints[previousI]);

        // Remove any small out-of-plane component before normalising.
        axis1 -= (normal & axis1)*normal;

        if (mag(axis1) < SMALL)
        {
            continue;
        }

        axis1 /= mag(axis1);
        const vector axis2(normal ^ axis1);

        scalar min1 = VGREAT;
        scalar max1 = -VGREAT;
        scalar min2 = VGREAT;
        scalar max2 = -VGREAT;

        forAll(facePoints, testPointI)
        {
            const vector d(facePoints[testPointI] - facePoints[previousI]);
            const scalar projection1 = d & axis1;
            const scalar projection2 = d & axis2;

            min1 = Foam::min(min1, projection1);
            max1 = Foam::max(max1, projection1);
            min2 = Foam::min(min2, projection2);
            max2 = Foam::max(max2, projection2);
        }

        const scalar area = (max1 - min1)*(max2 - min2);

        if (area < minArea)
        {
            minArea = area;
            bestAxis1 = axis1;
            bestAxis2 = axis2;
            bestCentre =
                facePoints[previousI]
              + 0.5*((min1 + max1)*axis1 + (min2 + max2)*axis2);
        }
    }

    if (minArea == VGREAT)
    {
        FatalErrorInFunction
            << "The face has no non-zero in-plane edge" << abort(FatalError);
    }

    axes_ = tensor(normal, bestAxis1, bestAxis2);
    centre_ =
        bestCentre + 0.5*(normalPosExtent - normalNegExtent)*normal;

    vector minLocal(VGREAT, VGREAT, VGREAT);
    vector maxLocal(-VGREAT, -VGREAT, -VGREAT);

    forAll(facePoints, pointI)
    {
        const point positivePoint
        (
            facePoints[pointI] + normalPosExtent*normal
        );
        const point negativePoint
        (
            facePoints[pointI] - normalNegExtent*normal
        );
        addToLocalBounds
        (
            axes_ & (positivePoint - centre_),
            minLocal,
            maxLocal
        );
        addToLocalBounds
        (
            axes_ & (negativePoint - centre_),
            minLocal,
            maxLocal
        );
    }

    halfLength_ = 0.5*scaleFactor*(maxLocal - minLocal);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OBB::OBB()
:
    centre_(point::zero),
    halfLength_(vector::zero),
    axes_(tensor::I)
{}


Foam::OBB::OBB(const point& centre, const vector& halfLength)
:
    centre_(centre),
    halfLength_(halfLength),
    axes_(tensor::I)
{}


Foam::OBB::OBB
(
    const point& centre,
    const vector& halfLength,
    const tensor& axes
)
:
    centre_(centre),
    halfLength_(halfLength),
    axes_(axes)
{
    orthonormaliseAxes();
}


Foam::OBB::OBB(const boundBox& box)
:
    centre_(box.midpoint()),
    halfLength_(0.5*(box.max() - box.min())),
    axes_(tensor::I)
{}


Foam::OBB::OBB(const pointField& points)
:
    centre_(point::zero),
    halfLength_(vector::zero),
    axes_(tensor::I)
{
    makeOBB(points);
}


Foam::OBB::OBB(const tmp<pointField>& points)
:
    centre_(point::zero),
    halfLength_(vector::zero),
    axes_(tensor::I)
{
    makeOBB(points());
    points.clear();
}


Foam::OBB::OBB
(
    const pointField& facePoints,
    const vector& faceNormal,
    const scalar normalPosExtent,
    const scalar normalNegExtent,
    const scalar scaleFactor
)
:
    centre_(point::zero),
    halfLength_(vector::zero),
    axes_(tensor::I)
{
    makeOBB
    (
        facePoints,
        faceNormal,
        normalPosExtent,
        normalNegExtent,
        scaleFactor
    );
}


Foam::OBB::OBB(Istream& is)
:
    centre_(point::zero),
    halfLength_(vector::zero),
    axes_(tensor::I)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OBB::overlaps(const OBB& box) const
{
    const vector& a = halfLength_;
    const vector& b = box.halfLength_;

    // Rotation from the other box coordinates into this box coordinates.
    const tensor relativeAxes(axes_ & box.axes_.T());
    tensor magRelativeAxes(tensor::zero);

    for (direction i = 0; i < vector::nComponents; ++i)
    {
        for (direction j = 0; j < vector::nComponents; ++j)
        {
            magRelativeAxes(i, j) = mag(relativeAxes(i, j));
        }
    }

    // Centre separation expressed in this box coordinates.
    const vector separation(axes_ & (box.centre_ - centre_));

    // The three axes of this box.
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        const scalar projectedB =
            b.x()*magRelativeAxes(i, 0)
          + b.y()*magRelativeAxes(i, 1)
          + b.z()*magRelativeAxes(i, 2);

        if (mag(separation[i]) > a[i] + projectedB)
        {
            return false;
        }
    }

    // The three axes of the other box.
    for (direction j = 0; j < vector::nComponents; ++j)
    {
        const scalar projectedSeparation = mag
        (
            separation.x()*relativeAxes(0, j)
          + separation.y()*relativeAxes(1, j)
          + separation.z()*relativeAxes(2, j)
        );
        const scalar projectedA =
            a.x()*magRelativeAxes(0, j)
          + a.y()*magRelativeAxes(1, j)
          + a.z()*magRelativeAxes(2, j);

        if (projectedSeparation > projectedA + b[j])
        {
            return false;
        }
    }

    // The nine cross products of one axis from each box.
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        const direction i1 = (i + 1) % vector::nComponents;
        const direction i2 = (i + 2) % vector::nComponents;

        for (direction j = 0; j < vector::nComponents; ++j)
        {
            // A parallel pair has no cross-product separating axis.
            if
            (
                1.0 - relativeAxes(i, j)*relativeAxes(i, j)
              <= parallelTolerance_*parallelTolerance_
            )
            {
                continue;
            }

            const direction j1 = (j + 1) % vector::nComponents;
            const direction j2 = (j + 2) % vector::nComponents;
            const scalar projectedSeparation = mag
            (
                separation[i2]*relativeAxes(i1, j)
              - separation[i1]*relativeAxes(i2, j)
            );
            const scalar projectedA =
                a[i1]*magRelativeAxes(i2, j)
              + a[i2]*magRelativeAxes(i1, j);
            const scalar projectedB =
                b[j1]*magRelativeAxes(i, j2)
              + b[j2]*magRelativeAxes(i, j1);

            if (projectedSeparation > projectedA + projectedB)
            {
                return false;
            }
        }
    }

    return true;
}


bool Foam::OBB::operator==(const OBB& box) const
{
    return
        centre_ == box.centre_
     && halfLength_ == box.halfLength_
     && axes_ == box.axes_;
}


// * * * * * * * * * * * * * * IOstream Operators * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const OBB& box)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << box.centre_ << token::SPACE
            << box.halfLength_ << token::SPACE
            << box.axes_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&box.centre_),
            sizeof(OBB)
        );
    }

    os.check("Ostream& operator<<(Ostream&, const OBB&)");
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, OBB& box)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> box.centre_ >> box.halfLength_ >> box.axes_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&box.centre_),
            sizeof(OBB)
        );
    }

    is.check("Istream& operator>>(Istream&, OBB&)");
    return is;
}


// ************************************************************************* //
