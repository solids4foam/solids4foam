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

#include "valveSliceAxis.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurface.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "IFstream.H"
#include "OSspecific.H"
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(valveSliceAxis, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        valveSliceAxis,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

Foam::fileName Foam::immersedBodyMotions::valveSliceAxis::surfaceFile
(
    const fileName& f
)
{
    fileName file(f);
    file.expand();

    if (!file.isAbsolute())
    {
        file = fileName("<constant>/triSurface"/file);
        file.expand();
    }

    if (!isFile(file))
    {
        FatalErrorInFunction
            << "Cannot find the file " << file << exit(FatalError);
    }

    return file;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBodyMotions::valveSliceAxis::calcFixed
(
    const pointField& points0
) const
{
    fixed_.setSize(points0.size(), false);
    fixedComputed_ = true;

    if (!annulusSurface_.empty())
    {
        const triSurface annulus(surfaceFile(annulusSurface_));
        const pointField& annPts = annulus.points();

        // Nearest annulus vertex to each valve vertex
        treeBoundBox bb(annPts);
        bb.inflate(1e-4);

        const indexedOctree<treeDataPoint> tree
        (
            treeDataPoint(annPts),
            bb,
            8,
            10,
            3.0
        );

        const scalar tol2 = sqr(annulusTolerance_);

        forAll(points0, i)
        {
            if (tree.findNearest(points0[i], tol2).hit())
            {
                fixed_[i] = true;
            }
        }

        const label nFixed = std::count(fixed_.begin(), fixed_.end(), true);

        Info<< "    valveSliceAxis: " << nFixed << " of " << points0.size()
            << " vertices within " << annulusTolerance_
            << " m of the annulus " << annulusSurface_ << " are fixed"
            << endl;

        if (nFixed == 0)
        {
            WarningInFunction
                << "No vertex of the valve surface is within "
                << annulusTolerance_ << " m of the annulus "
                << annulusSurface_ << endl;
        }
    }
    else if (!xiFile_.empty())
    {
        IFstream is(surfaceFile(xiFile_));

        // One value per vertex of the valve surface
        scalarField xi(points0.size());

        forAll(xi, i)
        {
            is >> xi[i];
        }

        token extra(is);
        if (extra.good())
        {
            FatalIOErrorInFunction(is)
                << xiFile_ << " has more values than the " << xi.size()
                << " vertices of the valve surface" << exit(FatalIOError);
        }

        // A field with the opposite sign convention
        label nNegative = 0;
        forAll(xi, i)
        {
            if (xi[i] < 0)
            {
                ++nNegative;
            }
        }

        const scalar sign =
            (max(xi) <= 0 && nNegative > xi.size()/2) ? -1 : 1;

        scalar minPositive = GREAT;
        scalar maxPositive = 0;

        forAll(xi, i)
        {
            xi[i] = max(sign*xi[i], 0);

            if (xi[i] > SMALL)
            {
                minPositive = min(minPositive, xi[i]);
                maxPositive = max(maxPositive, xi[i]);
            }
        }

        // As in the original: the values are scaled to a maximum of
        // activeLen, and a vertex is fixed where the scaled value is below
        // the smaller of 0.05 nearRamp and 0.25 nearRamp or half the smallest
        // positive value before scaling
        const scalar scale =
            (maxPositive > SMALL && activeLen_ > SMALL)
          ? activeLen_/maxPositive
          : 1;

        const scalar tol =
            min
            (
                max(SMALL, min(0.25*nearRamp_, 0.5*minPositive)),
                0.05*nearRamp_
            );

        forAll(xi, i)
        {
            fixed_[i] = (scale*xi[i] <= tol);
        }

        Info<< "    valveSliceAxis: "
            << std::count(fixed_.begin(), fixed_.end(), true) << " of "
            << points0.size() << " vertices are fixed by the values in "
            << xiFile_ << endl;
    }
}


Foam::scalar Foam::immersedBodyMotions::valveSliceAxis::spaceGain
(
    const scalar xif
) const
{
    if (xif < 0)
    {
        return 0;
    }
    else if (xif > nearRamp_)
    {
        return (xif <= activeLen_) ? 1 : 0;
    }

    const scalar s = xif/nearRamp_;

    if (spaceLaw_ == "linear")
    {
        return s;
    }
    else if (spaceLaw_ == "smoothstep")
    {
        return s*s*(3 - 2*s);
    }

    return s*s*s*(10 - 15*s + 6*s*s);
}


void Foam::immersedBodyMotions::valveSliceAxis::geometry
(
    const pointField& points0,
    vectorField& r0,
    scalarField& g
) const
{
    if (!fixedComputed_ || fixed_.size() != points0.size())
    {
        calcFixed(points0);
    }

    const vector axis(topCentre_ - bottomCentre_);
    const scalar L = mag(axis);
    const vector u(axis/L);

    r0.setSize(points0.size());
    g.setSize(points0.size());

    forAll(points0, i)
    {
        const scalar xi =
            min(max((points0[i] - bottomCentre_) & u, 0), L);

        r0[i] = points0[i] - (bottomCentre_ + xi*u);
        r0[i] -= (r0[i] & u)*u;

        g[i] = fixed_[i] ? 0 : spaceGain(fixedBottom_ ? xi : L - xi);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::valveSliceAxis::valveSliceAxis
(
    const dictionary& dict
)
:
    immersedBodyMotion(dict),
    timeLaw_(dict),
    bottomCentre_(dict.get<point>("bottomCentre")),
    topCentre_(dict.get<point>("topCentre")),
    fixedBottom_(true),
    Fmax_(dict.getCheck<scalar>("Fmax", scalarMinMax(0, 1))),
    spaceLaw_(dict.getOrDefault<word>("spaceLaw", "smootherstep")),
    nearRamp_(dict.getCheck<scalar>("nearRamp", scalarMinMax::ge(SMALL))),
    activeLen_(dict.get<scalar>("activeLen")),
    annulusSurface_(dict.getOrDefault<fileName>("annulusSurface", "")),
    annulusTolerance_
    (
        annulusSurface_.empty()
      ? 0
      : dict.getCheck<scalar>("annulusTolerance", scalarMinMax::ge(SMALL))
    ),
    xiFile_(dict.getOrDefault<fileName>("xiFile", "")),
    fixed_(),
    fixedComputed_(false)
{
    const word fixedEnd(dict.get<word>("fixedEnd"));

    if (fixedEnd == "top")
    {
        fixedBottom_ = false;
    }
    else if (fixedEnd != "bottom")
    {
        FatalIOErrorInFunction(dict)
            << "fixedEnd must be bottom or top, not " << fixedEnd
            << exit(FatalIOError);
    }

    if (mag(topCentre_ - bottomCentre_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "topCentre and bottomCentre must be different"
            << exit(FatalIOError);
    }

    if
    (
        spaceLaw_ != "linear"
     && spaceLaw_ != "smoothstep"
     && spaceLaw_ != "smootherstep"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown spaceLaw " << spaceLaw_ << nl
            << "Valid laws are: linear smoothstep smootherstep"
            << exit(FatalIOError);
    }

    if (!annulusSurface_.empty() && !xiFile_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Specify at most one of annulusSurface and xiFile"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::valveSliceAxis::points
(
    const pointField& points0,
    const scalar t
) const
{
    vectorField r0;
    scalarField g;
    geometry(points0, r0, g);

    const scalar G = timeLaw_.value(t);

    tmp<pointField> tpts(new pointField(points0));
    pointField& pts = tpts.ref();

    forAll(pts, i)
    {
        pts[i] -= min(Fmax_*g[i]*G, 1)*r0[i];
    }

    return tpts;
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::valveSliceAxis::velocity
(
    const pointField& x,
    const scalar t
) const
{
    FatalErrorInFunction
        << "The velocity of the deforming motion " << typeName
        << " is only defined at the surface points" << abort(FatalError);

    return tmp<vectorField>::New(x.size(), Zero);
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::valveSliceAxis::pointVelocities
(
    const pointField& points0,
    const scalar t
) const
{
    vectorField r0;
    scalarField g;
    geometry(points0, r0, g);

    const scalar G = timeLaw_.value(t);
    const scalar dGdt = timeLaw_.derivative(t);

    tmp<vectorField> tU(new vectorField(points0.size(), Zero));
    vectorField& U = tU.ref();

    forAll(U, i)
    {
        if (Fmax_*g[i]*G < 1)
        {
            U[i] = -Fmax_*g[i]*dGdt*r0[i];
        }
    }

    return tU;
}


// ************************************************************************* //
