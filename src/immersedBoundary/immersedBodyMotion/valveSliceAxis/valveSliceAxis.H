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
    Foam::immersedBodyMotions::valveSliceAxis

Description
    Deformation of the immersed surface of a heart valve by the radial
    collapse of its slices normal to an axis, from the reference (open)
    surface towards the axis (closed).

    The axis runs from \c bottomCentre \f$ C_b \f$ to \c topCentre
    \f$ C_t \f$, with the unit vector \f$ u = (C_t - C_b)/L \f$,
    \f$ L = |C_t - C_b| \f$. A surface point with reference position
    \f$ X_0 \f$ has the axial coordinate
    \f$ \xi = \min(\max((X_0 - C_b) \cdot u, 0), L) \f$, the centre of its
    slice \f$ C = C_b + \xi u \f$, and the radial vector \f$ r_0 \f$, the
    part of \f$ X_0 - C \f$ normal to the axis. Its position at time t is

    \f[
        X = X_0 - F r_0, \qquad
        F = \min(F_{max} \, g(\xi_f) \, G(t), 1)
    \f]

    where \f$ \xi_f \f$ is the distance along the axis from the fixed end
    (\c fixedEnd \c bottom, \f$ \xi_f = \xi \f$, or \c top,
    \f$ \xi_f = L - \xi \f$), G(t) the valve time law (see valveTimeLaw), and
    the space law g rises from 0 at the fixed end to 1 at \f$ \xi_f = \f$
    \c nearRamp, following \c spaceLaw (\c linear, \c smoothstep or
    \c smootherstep), is 1 up to \f$ \xi_f = \f$ \c activeLen, and 0 beyond.
    The velocity of the point is \f$ -F_{max} \, g(\xi_f) \, \dot{G}(t) \,
    r_0 \f$ while \f$ F < 1 \f$, and 0 once F is limited to 1.

    The points attached to the annulus can be held fixed, either those within
    \c annulusTolerance of a vertex of the surface \c annulusSurface, or
    those with a zero or small value in \c xiFile, a list of one scalar per
    vertex of the valve surface, in the order of its vertices, as in the
    original: the sign of the values is reversed if all are non-positive and
    most are negative, the negative values are set to zero, the values are
    scaled to a maximum of \c activeLen, and the vertices with a scaled value
    below the smaller of 0.05 \c nearRamp and 0.25 \c nearRamp or half the
    smallest positive value before scaling are fixed. Relative file names are
    relative to constant/triSurface.

    Usage:
    \verbatim
    motion
    {
        type            valveSliceAxis;
        bottomCentre    (0 0 0);        // [m]
        topCentre       (0 0 0.02);     // [m]
        fixedEnd        bottom;         // bottom or top
        Fmax            0.95;           // in [0, 1]
        spaceLaw        smootherstep;   // optional, default smootherstep
        nearRamp        0.004;          // [m]
        activeLen       0.016;          // [m]

        // Valve time law, see valveTimeLaw
        period          0.8;
        timeLaw         twoWindow;
        closeWindow     (0.625 0.655);
        openWindow      (0.94 1);

        // Optional: hold the vertices attached to the annulus fixed
        annulusSurface   "annulus.stl";
        annulusTolerance 1e-4;          // [m]
        // or
        xiFile           "valveXi.dat";
    }
    \endverbatim

    Without \c annulusSurface and \c xiFile this is the slice-axis valve
    motion of the immersed boundary library contributed to cardiacFoam by
    Sairam Pamulaparthi Venkata, and with one of them its saddle variant. The
    entries of the original are renamed (\c topCenter, \c bottomCenter,
    \c xiFieldFile, \c annulusStlFile and \c annulusTol become
    \c topCentre, \c bottomCentre, \c xiFile, \c annulusSurface and
    \c annulusTolerance), have no defaults, and the mode is set by the
    presence of \c annulusSurface or \c xiFile rather than by \c valveMode;
    the legacy axis from \c topCenter to \c bottomCenter with
    \c fixedEnd \c top is the same as the axis from \c bottomCentre to
    \c topCentre used here with \c fixedEnd \c top.

Author
    Based on the slice-axis valve motion of Sairam Pamulaparthi Venkata, UCD.
    Philip Cardiff, UCD.

SourceFiles
    valveSliceAxis.C

\*---------------------------------------------------------------------------*/

#ifndef valveSliceAxis_H
#define valveSliceAxis_H

#include "immersedBodyMotion.H"
#include "valveTimeLaw.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                       Class valveSliceAxis Declaration
\*---------------------------------------------------------------------------*/

class valveSliceAxis
:
    public immersedBodyMotion
{
    // Private Data

        //- Time law
        const valveTimeLaw timeLaw_;

        //- Ends of the axis [m]
        const point bottomCentre_;
        const point topCentre_;

        //- Is the fixed end the bottom?
        bool fixedBottom_;

        //- Maximum closure fraction
        const scalar Fmax_;

        //- Space law
        const word spaceLaw_;

        //- Length of the ramp of the space law from the fixed end [m]
        const scalar nearRamp_;

        //- Length of the active part from the fixed end [m]
        const scalar activeLen_;

        //- Surface of the annulus, optional
        const fileName annulusSurface_;

        //- Distance within which a point is on the annulus [m]
        const scalar annulusTolerance_;

        //- File of the axial coordinate of the vertices, optional
        const fileName xiFile_;

        //- Fixed vertices, computed on first use
        mutable boolList fixed_;

        //- Have the fixed vertices been computed?
        mutable bool fixedComputed_;


    // Private Member Functions

        //- Compute the fixed vertices of the reference surface points0
        void calcFixed(const pointField& points0) const;

        //- Space law at the distance xif from the fixed end
        scalar spaceGain(const scalar xif) const;

        //- Radial vectors and space gains of the reference points
        void geometry
        (
            const pointField& points0,
            vectorField& r0,
            scalarField& g
        ) const;


public:

    //- Runtime type information
    TypeName("valveSliceAxis");


    // Static Member Functions

        //- Expand the file name, relative to constant/triSurface
        static fileName surfaceFile(const fileName& f);


    // Constructors

        //- Construct from dictionary
        explicit valveSliceAxis(const dictionary& dict);


    //- Destructor
    virtual ~valveSliceAxis() = default;


    // Member Functions

        //- Positions at time t of the surface points
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const;

        //- The motion is not rigid
        virtual bool rigid() const
        {
            return false;
        }

        //- Not defined for a deforming motion: use pointVelocities
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const;

        //- Velocity at time t of the surface points whose reference
        //  positions are points0
        virtual tmp<vectorField> pointVelocities
        (
            const pointField& points0,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBodyMotions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
