
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFCoarsening_H
#define RBFCoarsening_H

#include "RBFInterpolation.H"
#include "fvCFD.H"
#include <fstream>

namespace rbf
{
    class RBFCoarsening
    {
        public:
            RBFCoarsening();

            explicit RBFCoarsening( std::shared_ptr<RBFInterpolation> rbf );

            RBFCoarsening(
                std::shared_ptr<RBFInterpolation> rbf,
                bool enabled,
                bool livePointSelection,
                bool livePointSelectionSumValues,
                scalar tol,
                scalar tolLivePointSelection,
                int coarseningMinPoints,
                int coarseningMaxPoints,
                bool exportTxt
                );

            RBFCoarsening(
                std::shared_ptr<RBFInterpolation> rbf,
                bool enabled,
                bool livePointSelection,
                bool livePointSelectionSumValues,
                scalar tol,
                scalar tolLivePointSelection,
                int coarseningMinPoints,
                int coarseningMaxPoints,
                bool twoPointSelection,
                bool exportTxt
                );

            RBFCoarsening(
                std::shared_ptr<RBFInterpolation> rbf,
                bool enabled,
                bool livePointSelection,
                bool livePointSelectionSumValues,
                scalar tol,
                scalar tolLivePointSelection,
                int coarseningMinPoints,
                int coarseningMaxPoints,
                bool twoPointSelection,
                bool surfaceCorrection,
                scalar ratioRadiusError,
                bool exportTxt
                );

            void greedySelection( const matrix & values );

            void compute(
                const matrix & positions,
                const matrix & positionsInterpolation
                );

            void interpolate(
                const matrix & values,
                matrix & valuesInterpolation
                );

            void setNbMovingAndStaticFaceCenters(
                int nbMovingFaceCenters,
                int nbStaticFaceCenters
                );

            void correctSurface( matrix & valuesInterpolation );

            std::shared_ptr<RBFInterpolation> rbf;
            std::shared_ptr<RBFInterpolation> rbfCoarse;
            bool enabled;
            bool livePointSelection;
            bool livePointSelectionSumValues;
            scalar tol;
            scalar tolLivePointSelection;
            int coarseningMinPoints;
            int coarseningMaxPoints;
            bool twoPointSelection;
            bool surfaceCorrection;
            scalar ratioRadiusError;
            bool exportTxt;
            Eigen::VectorXi selectedPositions;
            int nbStaticFaceCentersRemove;
            matrix positions;
            matrix positionsInterpolation;
            matrix values;
            matrix errorInterpolationCoarse;
            Eigen::VectorXi closestBoundaryIndexCorrection;
            matrix valuesCorrection;
            int nbMovingFaceCenters;
            int fileExportIndex;

            //static debug::debugSwitch debug;
    };
}

#endif
