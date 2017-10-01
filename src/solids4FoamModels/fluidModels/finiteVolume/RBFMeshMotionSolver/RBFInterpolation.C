
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"
#include "TPSFunction.H"

namespace rbf
{
    RBFInterpolation::RBFInterpolation()
        :
        rbfFunction( std::shared_ptr<RBFFunctionInterface> ( new TPSFunction() ) ),
        polynomialTerm( true ),
        cpu( false ),
        computed( false ),
        n_A( 0 ),
        n_B( 0 ),
        dimGrid( 0 ),
        Hhat(),
        Phi(),
        lu(),
        positions(),
        positionsInterpolation()
    {}

    RBFInterpolation::RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction )
        :
        rbfFunction( rbfFunction ),
        polynomialTerm( true ),
        cpu( false ),
        computed( false ),
        n_A( 0 ),
        n_B( 0 ),
        dimGrid( 0 ),
        Hhat(),
        Phi(),
        lu(),
        positions(),
        positionsInterpolation()
    {
        assert( rbfFunction );
    }

    RBFInterpolation::RBFInterpolation(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        bool polynomialTerm,
        bool cpu
        )
        :
        rbfFunction( rbfFunction ),
        polynomialTerm( polynomialTerm ),
        cpu( cpu ),
        computed( false ),
        n_A( 0 ),
        n_B( 0 ),
        dimGrid( 0 ),
        Hhat(),
        Phi(),
        lu(),
        positions(),
        positionsInterpolation()
    {
        assert( rbfFunction );
    }

    void RBFInterpolation::evaluateH(
        const matrix & positions,
        matrix & H
        )
    {
        // RBF function evaluation

        scalar r = 0;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = i; j < n_A; j++ )
            {
                r = ( positions.row( i ) - positions.row( j ) ).norm();
                H( j, i ) = rbfFunction->evaluate( r );
            }
        }
    }

    void RBFInterpolation::evaluatePhi(
        const matrix & positions,
        const matrix & positionsInterpolation,
        matrix & Phi
        )
    {
        // Evaluate Phi which contains the evaluation of the radial basis function

        scalar r = 0;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_B; j++ )
            {
                r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
                Phi( j, i ) = rbfFunction->evaluate( r );
            }
        }
    }

    void RBFInterpolation::compute(
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        // Verify input

        assert( positions.cols() == positionsInterpolation.cols() );
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( positionsInterpolation.rows() > 0 );

        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        dimGrid = positions.cols();

        // Radial basis function interpolation
        // Initialize matrices H and Phi
        matrix H( n_A, n_A ), Phi( n_B, n_A );

        if ( polynomialTerm )
        {
            H.resize( n_A + dimGrid + 1, n_A + dimGrid + 1 );
        }

        // Evaluate radial basis functions for matrix H

        evaluateH( positions, H );

        // Include polynomial contributions
        if ( polynomialTerm )
        {
            for ( int i = 0; i < n_A; i++ )
                H( n_A, i ) = 1;

            H.bottomLeftCorner( dimGrid, n_A ) = positions.block( 0, 0, n_A, dimGrid ).transpose();

            for ( int i = 0; i < dimGrid + 1; i++ )
                for ( int j = 0; j < dimGrid + 1; j++ )
                    H( H.rows() - dimGrid - 1 + i, H.rows() - dimGrid - 1 + j ) = 0;
        }

        if ( cpu )
        {
            this->positions = positions;
            this->positionsInterpolation = positionsInterpolation;

            lu.compute( H.selfadjointView<Eigen::Lower>() );
        }

        if ( not cpu )
        {
            // Evaluate Phi which contains the evaluation of the radial basis function

            if ( polynomialTerm )
                Phi.resize( n_B, n_A + dimGrid + 1 );

            evaluatePhi( positions, positionsInterpolation, Phi );

            // Include polynomial contributions in matrix Phi
            if ( polynomialTerm )
            {
                for ( int i = 0; i < Phi.rows(); i++ )
                    Phi( i, n_A ) = 1;

                Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );
            }

            // Compute the LU decomposition of the matrix H
            Eigen::FullPivLU<matrix> lu( H.selfadjointView<Eigen::Lower>() );

            // Compute interpolation matrix

            Hhat.noalias() = Phi * lu.inverse();

            Hhat.conservativeResize( n_B, n_A );
        }

        computed = true;
    }

    void RBFInterpolation::interpolate(
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        if ( cpu && not computed )
            compute( positions, positionsInterpolation );

        assert( computed );

        if ( cpu )
        {
            matrix B, valuesLU( n_A, values.cols() ), Phi( n_B, n_A );

            if ( polynomialTerm )
            {
                Phi.resize( n_B, n_A + dimGrid + 1 );
                valuesLU.resize( n_A + dimGrid + 1, values.cols() );
            }

            valuesLU.setZero();
            valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;

            B = lu.solve( valuesLU );

            evaluatePhi( positions, positionsInterpolation, Phi );

            if ( polynomialTerm )
            {
                // Include polynomial contributions in matrix Phi

                for ( int i = 0; i < Phi.rows(); i++ )
                    Phi( i, n_A ) = 1;

                Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );
            }

            valuesInterpolation.noalias() = Phi * B;
        }

        if ( not cpu )
        {
            valuesInterpolation.noalias() = Hhat * values;
        }

        assert( valuesInterpolation.rows() == n_B );
        assert( values.cols() == valuesInterpolation.cols() );
    }

    /*
     * Compute interpolation matrix and directly interpolate the values.
     * The algorithms solves for the coefficients, and explicitly
     * uses the coefficients to interpolate the data to the new positions.
     */
    void RBFInterpolation::interpolate(
        const matrix & positions,
        const matrix & positionsInterpolation,
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        // Verify input

        assert( positions.cols() == positionsInterpolation.cols() );
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( positionsInterpolation.rows() > 0 );

        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        dimGrid = positions.cols();

        // Radial basis function interpolation

        // Initialize matrices

        matrix H( n_A, n_A );

        if ( polynomialTerm )
            H.resize( n_A + dimGrid + 1, n_A + dimGrid + 1 );

        // RBF function evaluation

        evaluateH( positions, H );

        // Include polynomial contributions

        if ( polynomialTerm )
        {
            // THIJS: initialize Phi if empty
            if ( Phi.cols() == 0 )
            {
                Phi.conservativeResize( n_B, dimGrid + 1 );
            }

            for ( int i = 0; i < n_A; i++ )
                H( n_A, i ) = 1;

            H.bottomLeftCorner( dimGrid, n_A ) = positions.block( 0, 0, n_A, dimGrid ).transpose();

            for ( int i = 0; i < dimGrid + 1; i++ )
                for ( int j = 0; j < dimGrid + 1; j++ )
                    H( H.rows() - dimGrid - 1 + i, H.rows() - dimGrid - 1 + j ) = 0;
        }

        // Calculate coefficients gamma and beta

        matrix B;

        matrix valuesLU( H.rows(), values.cols() );
        valuesLU.setZero();
        valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;

        lu.compute( H.selfadjointView<Eigen::Lower>() );
        B = lu.solve( valuesLU );

        // Evaluate Phi_BA which contains the evaluation of the radial basis function
        // This method is only used by the greedy algorithm, and the matrix Phi
        // is therefore enlarged at every greedy step.

        buildPhi( positions, positionsInterpolation );

        if ( polynomialTerm )
        {
            // Include polynomial contributions in matrix Phi

            for ( int i = 0; i < Phi.rows(); i++ )
                Phi( i, n_A ) = 1;

            Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );
        }

        valuesInterpolation.noalias() = Phi * B;

        computed = true;
    }

    void RBFInterpolation::buildPhi(
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        int phiColsOld = Phi.cols();

        if ( polynomialTerm )
            Phi.conservativeResize( n_B, n_A + dimGrid + 1 );
        else
            Phi.conservativeResize( n_B, n_A );

        int nNewPoints = Phi.cols() - phiColsOld;

        if ( nNewPoints == Phi.cols() )
            nNewPoints = n_A;

        scalar r = 0;

        for ( int i = 0; i < nNewPoints; i++ )
        {
            int index = Phi.cols() - (i + 1);

            if ( polynomialTerm )
                index = Phi.cols() - 1 - dimGrid - (i + 1);

            for ( int j = 0; j < n_B; j++ )
            {
                r = ( positions.row( index ) - positionsInterpolation.row( j ) ).norm();
                Phi( j, index ) = rbfFunction->evaluate( r );
            }
        }
    }

    /*
     * This function is only called by the RBFCoarsening class.
     * It is assumed that the polynomial term is included in the
     * interpolation, and that the fullPivLu decomposition is
     * used to solve for the coefficients B.
     */
    void RBFInterpolation::interpolate2(
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        assert( computed );

        matrix valuesLU( values.rows(), values.cols() );

        // resize valuesLU if polynomial is used
        if ( polynomialTerm )
        {
            valuesLU.conservativeResize( values.rows() + values.cols() + 1, values.cols() );
        }

        valuesLU.setZero(); // initialize all values zero

        // Set correct part of valuesLU equal to values
        if ( polynomialTerm )
        {
            valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;
        }
        else
        {
            valuesLU = values;
        }

        valuesInterpolation.noalias() = Phi * lu.solve( valuesLU );

        assert( valuesInterpolation.rows() == n_B );
        assert( values.cols() == valuesInterpolation.cols() );
    }
}
