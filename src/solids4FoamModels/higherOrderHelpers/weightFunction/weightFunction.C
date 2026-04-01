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

#include "weightFunction.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#ifdef OPENFOAM_COM
const Enum<weightFunction::weightFunctionType>
weightFunction::weightFunctionNames_
({
    {weightFunction::weightFunctionType::ONE, "one"},
    {weightFunction::weightFunctionType::LINEAR, "linear"},
    {weightFunction::weightFunctionType::INV_DIST, "inverseDistance"},
    {
        weightFunction::weightFunctionType::RAD_SYMM_EXP,
        "radiallySymmetricExponential"
    }
});
#else
template<>
const char* NamedEnum<weightFunction::weightFunctionType, 4>::names[] =
{
    "one",
    "linear",
    "inverseDistance",
    "radiallySymmetricExponential"
};

const NamedEnum<weightFunction::weightFunctionType, 4>
weightFunction::weightFunctionNames_;
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

weightFunction::weightFunction
(
    const dictionary& dict
)
:
    type_(weightFunctionNames_[dict.lookupOrDefault<word>("type", "one")]),
    k_(dict.lookupOrDefault<scalar>("k", 6.0)),
    expSqrK_(Foam::exp(-sqr(k_))),
    invDenom_(1.0/(1.0 - expSqrK_)),
    s_(dict.lookupOrDefault<scalar>("s", 1000)),
    b_(dict.lookupOrDefault<scalar>("b", 3))
{
    if (type() == weightFunctionType::INV_DIST && (s_ <= 0 || b_ <= 0))
    {
        FatalErrorInFunction
            << "s and b parameters should be positive for "
            << weightFunctionNames_[type()] << ", but s is set to " << s_
            << " and b to " << b_
            << abort(FatalError);
    }

    if (type() == weightFunctionType::RAD_SYMM_EXP && k_ <= 0)
    {
        FatalErrorInFunction
            << "Kernel parameter k  must be positive for "
            << weightFunctionNames_[type()] << ", but k is set to " << k_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar weightFunction::weight(const scalar d, const scalar maxDist) const
{
    scalar w = -1;

    if (type_ == weightFunctionType::ONE)
    {
        w = 1.0;
    }
    else if (type_ == weightFunctionType::LINEAR)
    {
        w = 1 - (d/maxDist);
    }
    else if (type_ == weightFunctionType::INV_DIST)
    {
        w = 1.0/(1.0 + s_*pow((d/(2*maxDist)),b_));
    }
    else if (type_ == weightFunctionType::RAD_SYMM_EXP)
    {
        // Smoothing length
        const scalar dm = 2*maxDist;
        const scalar sqrK = -sqr(k_);

        w = (Foam::exp(sqr(d/dm)*sqrK) - expSqrK_)*invDenom_;
    }
    else
    {
        FatalErrorInFunction
            << "Unrecognised weight function. Available options are: "
            << weightFunctionNames_.sortedToc()
            << abort(FatalError);
    }

    // Clip small negative value
    w = max(SMALL, w);

    return w;
}

} // End namespace Foam

// ************************************************************************* //
