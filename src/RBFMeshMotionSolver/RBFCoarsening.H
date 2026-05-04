
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
    RBFCoarsening

Description
    Radial basis function coarsening utilities.

Author
    David Blom, TU Delft. All rights reserved.
    NOTE: This file is distributed under the GNU General Public License (GPL).
    The phrase "All rights reserved" above is retained from the original
    source but does not limit the permissions granted by the GPL.

SourceFiles
    RBFCoarsening.C

\*---------------------------------------------------------------------------*/

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
