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

Description
    Function that creates and writes principal stress fields given a stress
    tensor.

    Six fields are created:
      - sigmaMin, sigmaMid and sigmaMax: principal stresses ordered from most
        negative to most positive
      - sigmaMinDir, sigmaMidDir and sigmaMaxDir: corresponding principal
        directions

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef principalStressFields_H
#define principalStressFields_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Set principal stress values and directions
    void calculateEigenValues
    (
        const symmTensor sigma,
        const bool& compressionPositive,
        scalar& sigmaMax,
        scalar& sigmaMid,
        scalar& sigmaMin,
        vector& sigmaMaxDir,
        vector& sigmaMidDir,
        vector& sigmaMinDir
    );

    // Calculate and write principal stress value and direction fields
    void writePrincipalStressFields
    (
        const volSymmTensorField& sigma, const bool& compressionPositive
    );

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
