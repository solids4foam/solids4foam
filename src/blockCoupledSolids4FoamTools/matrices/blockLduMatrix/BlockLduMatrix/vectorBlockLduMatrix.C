/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#ifndef vectorBlockLduMatrix_H
#define vectorBlockLduMatrix_H

#include "BlockLduMatrix.H"
#include "solidPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// template<>
// void Foam::BlockLduMatrix<Foam::vector>::sumDiag()
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     this->decoupledSumDiag();
// }


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::negSumDiag()
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     this->decoupledNegSumDiag();
// }


// template<>
// void Foam::BlockLduMatrix<vector>::check() const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     this->decoupledCheck();
// }


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::relax
// (
//     const vectorField& x,
//     vectorField& b,
//     const scalar alpha
// )
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     this->decoupledRelax(x, b, alpha);
// }


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::operator*=(const scalarField& sf)
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     this->decoupledMultEqOp(sf);
// }


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::AmulCore
// (
//     vectorField& Ax,
//     const vectorField& x
// ) const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     //decoupledAmulCore(Ax, x);
// }


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::TmulCore
// (
//     vectorField& Tx,
//     const vectorField& x
// ) const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     //decoupledTmulCore(Tx, x);
// }


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::segregateB
// (
//     vectorField&,
//     const vectorField&
// ) const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     FatalErrorIn
//     (
//         "void Foam::BlockLduMatrix<vector>::segregateB\n"
//         "(\n"
//         "    vectorField&,\n"
//         "    const vectorField&\n"
//         ") const"
//     )   << "Requested decoupling of vector matrix - never coupled"
//         << abort(FatalError);
// }


// template<>
// Foam::tmp<Foam::vectorField>
// Foam::BlockLduMatrix<Foam::vector>::H(const vectorField& x) const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     return decoupledH(x);
// }


// template<>
// Foam::tmp<Foam::vectorField>
// Foam::BlockLduMatrix<Foam::vector>::faceH(const vectorField& x) const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     // Decoupled version
//     return decoupledFaceH(x);
// }

template<>
void Foam::BlockLduMatrix<Foam::vector>::Amul
(
    Field<vector>& Ax,
    const Field<vector>& x
) const
{
    Ax = pTraits<vector>::zero;

    // Initialise the update of coupled interfaces
    initInterfaces(coupleUpper_, Ax, x);

    AmulCore(Ax, x);

    // Update coupled interfaces
    updateInterfaces(coupleUpper_, Ax, x);

    // Philipc: Update global fields
    // This could be rewritten as a global patch, so no change would be required
    // here: but I get MPI errors when I call updateGlobalFields within the
    // global patch...
    if (Pstream::parRun())
    {
        mesh().thisDb().lookupObject<solidPolyMesh>
        (
            "solidPolyMesh"
        ).updateGlobalFields(Ax, x);
    }
}


// template<>
// void Foam::BlockLduMatrix<Foam::vector>::Tmul
// (
//     Field<vector>& Ax,
//     const Field<vector>& x
// ) const
// {
//     FatalErrorIn("myfunc")
//         << abort(FatalError);

//     FatalErrorIn
//     (
//         "void BlockLduMatrix<vector>::Tmul"
//     ) << "philipc: Tmul not implemented" << abort(FatalError);
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
