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
    adaptiveStencil

Description
    Helper class for the adaptive enlargement of least-squares stencils.

    All stencils are first built with the size requested by the
    reconstruction scheme. The condition number of each least-squares system
    is then estimated and the stencils with an estimate larger than
    maxConditionNumber are enlarged by extraCellsPerPass cells. This is
    repeated for the enlarged stencils until the estimate is acceptable or
    the stencil has grown by maxExtraCells cells.

    The number of cells requested for a stencil never exceeds maxExtraCells.
    Cells at the same distance as the last requested cell are always added
    together, which keeps the stencils symmetric on regular meshes, so a
    stencil can exceed maxExtraCells by the remaining cells of that last
    group of equally distant cells.

    The condition number is estimated as the ratio of the largest and the
    smallest diagonal entry of R from the column-pivoted QR decomposition of
    the weighted least-squares matrix. This is a lower bound of the
    2-norm condition number, typically within a factor of two, and it is much
    cheaper than the singular value decomposition. By default each polynomial
    term of the least-squares matrix is scaled to unit norm first, which
    removes the contribution of the Taylor series scaling to the condition
    number; the estimate then measures the quality of the stencil geometry
    only and is smaller than the condition number of the unscaled matrix.

    The class is activated by the optional adaptiveStencilCoeffs dictionary:

    \verbatim
        adaptiveStencilCoeffs
        {
            // Stencils with a larger estimate are enlarged
            maxConditionNumber  50;

            // Cells added to a stencil in one pass
            // Default: half of the number of polynomial terms
            extraCellsPerPass   4;

            // Stencils stop growing after this many cells are added
            // Default: three times the number of polynomial terms
            maxExtraCells       30;

            // Scale the polynomial terms before estimating
            // Default: yes
            normalise           yes;
        }
    \endverbatim

SourceFiles
    adaptiveStencil.C
    adaptiveStencilTemplates.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef adaptiveStencil_H
#define adaptiveStencil_H

#include "dictionary.H"
#include "Switch.H"
#include "boolList.H"
#include "labelList.H"
#include "scalarField.H"
#include "word.H"
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class adaptiveStencil Declaration
\*---------------------------------------------------------------------------*/

class adaptiveStencil
{
    // Private data

        //- Is the adaptive enlargement active
        bool active_;

        //- Stencils with a larger condition number estimate are enlarged
        scalar maxConditionNumber_;

        //- Number of cells added to a stencil in one pass
        label extraCellsPerPass_;

        //- Stencils stop growing after this many cells are added
        label maxExtraCells_;

        //- Scale the polynomial terms to unit norm before estimating
        Switch normalise_;

public:

    // Public data types

        //- Request for condition number estimates, passed to the function
        //  that assembles the least-squares systems. Estimates are
        //  calculated for the active stencils only.
        struct conditioning
        {
            //- Stencils for which the estimate is required
            const boolList& active;

            //- Condition number estimates
            scalarField& estimate;

            conditioning(const boolList& a, scalarField& e)
            :
                active(a),
                estimate(e)
            {}
        };


    // Constructors

        //- Construct from the reconstruction scheme dictionary and the number
        //  of polynomial terms
        adaptiveStencil(const dictionary& schemeDict, const label nTerms);


    // Member Functions

        //- Is the adaptive enlargement active
        bool active() const
        {
            return active_;
        }

        //- Stencils with a larger condition number estimate are enlarged
        scalar maxConditionNumber() const
        {
            return maxConditionNumber_;
        }

        //- Estimate the condition number of the weighted least-squares
        //  matrix. Polynomial terms are stored in the rows of Q and stencil
        //  points in the columns
        scalar conditionEstimate
        (
            const Eigen::MatrixXd& Q,
            const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& W
        ) const;

        //- Enlarge the stencils with a poor condition number estimate.
        //  Returns the estimates for the final stencils.
        //    calcEstimate(conditioning&) calculates the estimates of the
        //        active stencils
        //    stencilSize(label) returns the current size of a stencil
        //    enlarge(const labelUList&, const labelUList&) enlarges the
        //        given stencils by the given number of cells
        template<class CalcEstimate, class StencilSize, class Enlarge>
        scalarField adapt
        (
            const word& stencilType,
            const label nStencils,
            const CalcEstimate& calcEstimate,
            const StencilSize& stencilSize,
            const Enlarge& enlarge
        ) const;

        //- Write a field defined on cells
        static void writeCellField
        (
            const fvMesh& mesh,
            const word& name,
            const scalarField& values
        );

        //- Write a field defined on all mesh faces
        static void writeFaceField
        (
            const fvMesh& mesh,
            const word& name,
            const scalarField& values
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adaptiveStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
