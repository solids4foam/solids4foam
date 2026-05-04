
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
    RBFInterpolation

Description
    Radial basis function interpolation utilities.

Author
    David Blom, TU Delft. All rights reserved.
    NOTE: This file is distributed under the GNU General Public License (GPL).
    The phrase "All rights reserved" above is retained from the original
    source but does not limit the permissions granted by the GPL.

SourceFiles
    RBFInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef RBFInterpolation_H
#define RBFInterpolation_H

#include <memory>
#include <Eigen/Dense>
#include "RBFFunctionInterface.H"
#include "fvCFD.H"

namespace rbf
{
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> matrix;
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector;

    class RBFInterpolation
    {
        public:
            RBFInterpolation();

            explicit RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction );

            RBFInterpolation(
                std::shared_ptr<RBFFunctionInterface> rbfFunction,
                bool polynomialTerm,
                bool cpu
                );

            void compute(
                const matrix & positions,
                const matrix & positionsInterpolation
                );

            void interpolate(
                const matrix & values,
                matrix & valuesInterpolation
                );

            void interpolate(
                const matrix & positions,
                const matrix & positionsInterpolation,
                const matrix & values,
                matrix & valuesInterpolation
                );

            void interpolate2(
                const matrix & values,
                matrix & valuesInterpolation
                );

            void buildPhi(
                const matrix & positions,
                const matrix & positionsInterpolation
                );

            std::shared_ptr<RBFFunctionInterface> rbfFunction;
            bool polynomialTerm;
            bool cpu;
            bool computed;
            int n_A;
            int n_B;
            int dimGrid;
            matrix Hhat;
            matrix Phi;
            Eigen::FullPivLU<matrix> lu;
            matrix positions;
            matrix positionsInterpolation;

        private:
            void evaluateH(
                const matrix & positions,
                matrix & H
                );

            void evaluatePhi(
                const matrix & positions,
                const matrix & positionsInterpolation,
                matrix & Phi
                );
    };
}

#endif
