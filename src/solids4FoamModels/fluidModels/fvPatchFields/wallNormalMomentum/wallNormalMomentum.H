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

Namespace
    Foam::wallNormalMomentum

Description
    The viscous and source terms of the wall-normal momentum equation,
    shared by the wall pressure conditions movingWallPressure and
    elasticWallPressure.

    At a no-slip wall the convective term vanishes, so the normal momentum
    equation gives

        n & grad(p) = n & (viscous term + source - a)

    with the viscous term of the incompressible momentum equation

        laplacian(nuEff, U) + div(nuEff*dev2(T(grad(U))))

    and, on OpenFOAM.com, the momentum source of the fvOptions (e.g. a body
    force), evaluated in the cells next to the wall at the current iterate,
    with nuEff from the fluid's turbulence (momentum transport) model. The
    transpose part is the Gauss linear divergence of linearly interpolated
    face values, so that it needs no scheme of its own. The fvModels of
    OpenFOAM.org are not included.

    The evaluation uses finite volume operators on the whole field, i.e.
    processor-boundary communication: all processors must call it together.

SourceFiles
    wallNormalMomentum.C

\*---------------------------------------------------------------------------*/

#ifndef wallNormalMomentum_H
#define wallNormalMomentum_H

#include "volFields.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wallNormalMomentum
{

//- Return the viscous and source terms of the momentum equation of the
//  velocity field UName in all cells; patchName is used in error messages
tmp<volVectorField> field
(
    const objectRegistry& db,
    const word& UName,
    const word& patchName
);

//- Return the normal component of the given term in the cells next to the
//  patch
tmp<scalarField> patchNormal
(
    const fvPatch& patch,
    const volVectorField& term
);

} // End namespace wallNormalMomentum

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
