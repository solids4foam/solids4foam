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

#ifdef OPENFOAM_COM

#include "realEigenValues.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

vector realEigenValues(const tensor& t)
{
    const Vector<complex> compEig = eigenValues(t);

    return vector
    (
        compEig[vector::X].Re(),
        compEig[vector::Y].Re(),
        compEig[vector::Z].Re()
    );
}


tensor realEigenVectors(const tensor& t)
{
    const Tensor<complex> comEigVecs = eigenVectors(t);

    return tensor
    (
        comEigVecs[tensor::XX].Re(),
        comEigVecs[tensor::XY].Re(),
        comEigVecs[tensor::XZ].Re(),
        comEigVecs[tensor::YX].Re(),
        comEigVecs[tensor::YY].Re(),
        comEigVecs[tensor::YZ].Re(),
        comEigVecs[tensor::ZX].Re(),
        comEigVecs[tensor::ZY].Re(),
        comEigVecs[tensor::ZZ].Re()
    );
}

} // End namespace Foam

#endif // OPENFOAM_COM

// ************************************************************************* //
