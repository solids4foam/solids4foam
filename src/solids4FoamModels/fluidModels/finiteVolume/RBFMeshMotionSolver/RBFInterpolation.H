
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

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
