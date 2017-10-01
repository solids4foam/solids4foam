/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidPolyBoundaryMesh.H"
#include "polyBoundaryMesh.H"
//#include "faceTetPolyPatch.H"
//#include "globalTetPolyPatch.H"

#include "processorPolyPatch.H"

// #include "solidPolyMesh.H"
// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void solidPolyBoundaryMesh::addPatches(const polyBoundaryMesh& basicBdry)
// {
//     setSize(basicBdry.size());

//     // Set boundary patches
//     fvPatchList& Patches = *this;

//     forAll(Patches, patchI)
//     {
//         Patches.set(patchI, fvPatch::New(basicBdry[patchI], *this));
//     }
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyBoundaryMesh
solidPolyBoundaryMesh::solidPolyBoundaryMesh
(
    const solidPolyMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    //fvPatchList(basicBdry.size()),
    mesh_(m),
    polyBdry_(basicBdry)
{
    //addPatches(basicBdry);

    // if (Pstream::parRun())
    // {
    //     //FatalErrorIn
    //     WarningIn
    //     (
    //         "solidPolyBoundaryMesh::solidPolyBoundaryMesh"
    //     )   << " lduInterface not implemented yet for coupled solid solver"
    //         //<< abort(FatalError);
    //         << endl;
    // }

    // PC: Do we need to create special polyPatches...? or can we just use
    // polyBoundaryMesh functionality...

    // Set boundary patches
    //solidPolyPatchList& Patches = *this;

    // forAll(Patches, patchI)
    // {
    //     Patches.set
    //     (
    //         patchI,
    //         faceTetPolyPatch::New(basicBdry[patchI], *this).ptr()
    //     );
    // }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduInterfacePtrsList solidPolyBoundaryMesh::interfaces() const
{
    //lduInterfacePtrsList interfaces(size());
    //lduInterfacePtrsList interfaces(0);
    lduInterfacePtrsList interfaces(polyBdry_.size());

    // For now, we will just use the interfaces from the fvMesh
    // This is not the nicest way to do this but OK for now
    // const fvMesh& fMesh = mesh_.time().lookupObject<fvMesh>("region0");
    // const fvBoundaryMesh& fBdryMesh = fMesh.boundary();

    // forAll(interfaces, patchi)
    // {
    //     if (isA<lduInterface>(fBdryMesh[patchi]))
    //     {
    //         interfaces.set
    //         (
    //             patchi,
    //            &refCast<const lduInterface>(fBdryMesh[patchi])
    //         );
    //     }
    // }

    // Note: interfaces is overwritten with D.boundaryField().blockInterfaces()
    // during the discretisation
    return interfaces;
}


// const globalTetPolyPatch&
// solidPolyBoundaryMesh::globalPatch() const
// {
//     const solidPolyPatchList& patches = *this;

//     forAll (patches, patchI)
//     {
//         if (isA<globalTetPolyPatch>(patches[patchI]))
//         {
//             return refCast<const globalTetPolyPatch>
//             (
//                 patches[patchI]
//             );
//         }
//     }

//     FatalErrorIn
//     (
//         "const globalTetPolyPatch&"
//         "solidPolyBoundaryMesh::globalPatch() const"
//     )   << "patch not found.  Is this case running in parallel?"
//         << abort(FatalError);

//     // Dummy return
//     return refCast<const globalTetPolyPatch>(patches[0]);
// }


// faceListList solidPolyBoundaryMesh::boundaryTriFaces() const
// {
//     faceListList result(size());

//     forAll (result, patchI)
//     {
//         result[patchI] = operator[](patchI).triFaces();
//     }

//     return result;
// }


// void solidPolyBoundaryMesh::updateMesh()
// {
//     solidPolyPatchList& Patches = *this;

//     forAll(Patches, patchI)
//     {
//         Patches[patchI].updateMesh();
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
