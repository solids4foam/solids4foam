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

#include "gradientSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

gradientSchemes::gradientSchemes
(
    const fvMesh& vm
)
:
    mesh_(vm),
    own_(mesh_.owner()),
    nei_(mesh_.neighbour()),
    X_(mesh_.C()),
    XF_(mesh_.Cf()),
    XN_(mesh_.points()),

    Ainv_
    (
        IOobject("Ainv", mesh_),
        mesh_,
        dimensionedTensor("Ainv", dimensionSet(0,2,0,0,0,0,0), tensor::zero)
    ),

    AinvLocal_
    (
        IOobject("AinvLocal", mesh_),
        mesh_,
        dimensionedTensor
        (
            "AinvLocal",
            dimensionSet(0,2,0,0,0,0,0),
            tensor::zero
        )
    )
{
    gradientSchemes::distanceMatrix(Ainv_);
    gradientSchemes::distanceMatrixLocal(AinvLocal_);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

gradientSchemes::~gradientSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrix
(
    GeometricField<tensor, fvPatchField, volMesh>& U
)
{
    forAll(own_, faceID)
    {
        const label& ownCellID = own_[faceID];
        const label& neiCellID = nei_[faceID];
        const vector& dOwn = X_[neiCellID] - X_[ownCellID];
        const vector& dNei  = X_[ownCellID] - X_[neiCellID];

        U[ownCellID] += dOwn*dOwn;
        U[neiCellID] += dNei*dNei;
    }

    if (Pstream::parRun())
    {
        forAll(mesh_.boundary(), patchID)
        {
            if (mesh_.boundary()[patchID].coupled())
            {
                const vectorField X_nei
                (
                  X_.boundaryField()[patchID].patchNeighbourField()
                );

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = X_nei[facei] - X_[bCellID];
                    U[bCellID] += d*d;
                }
            }
        }

        U.correctBoundaryConditions();
    }

    U.primitiveFieldRef() = inv(U.internalField());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrixLocal
(
    GeometricField<tensor, fvPatchField, volMesh>& Ainv
) const
{
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tvf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject("distanceMatrixLocal", mesh_),
            mesh_,
            dimensioned<tensor>("0", Ainv.dimensions(), pTraits<tensor>::zero)
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> dCd = tvf();

    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];
        const vector& dOwn = XF_[faceID] - X_[ownID];
        const vector& dNei = XF_[faceID] - X_[neiID];

        dCd[ownID] += dOwn*dOwn;
        dCd[neiID] += dNei*dNei;
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = XF_.boundaryField()[patchID][facei] - X_[bCellID];
            dCd[bCellID] += d*d;

            // if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            // {
            //     const label& faceID =
            //         mesh_.boundary()[patchID].start() + facei;

            //     forAll(mesh_.faces()[faceID], nodei)
            //     {
            //         const label& nodeID = mesh_.faces()[faceID][nodei];

            //         d = XN_[nodeID] - X_[bCellID];
            //         dCd[bCellID] += d * d;

            //         for (int i=0; i<7; i++)
            //         {
            //             d =
            //                 ((((i+1)*XN_[nodeID])
            //               + ((7 - i)*XF_.boundaryField()[patchID][facei]))/8.0)
            //               - X_[bCellID];
            //             dCd[bCellID] += d * d;
            //         }
            //     }
            // }
        }
    }

    Ainv.primitiveFieldRef() = inv(dCd.internalField());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volVectorField gradientSchemes::gradient
(
    const GeometricField<scalar, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ugrad = tvf();

    forAll(own_, faceID)
    {
        const label& cellID = own_[faceID];
        const label& neiID = nei_[faceID];
        const vector& dcell = X_[neiID] - X_[cellID];
        const vector& dnei = X_[cellID] - X_[neiID];

        Ugrad[cellID] += Ainv_[cellID] & (U[neiID] - U[cellID])*dcell;
        Ugrad[neiID] += Ainv_[neiID] & (U[cellID] - U[neiID])*dnei;
    }

    if (Pstream::parRun())
    {
        forAll(mesh_.boundary(), patchID)
        {
            if (mesh_.boundary()[patchID].coupled())
            {
                const vectorField X_nei
                (
                  X_.boundaryField()[patchID].patchNeighbourField()
                );

                const scalarField U_nei
                (
                  U.boundaryField()[patchID].patchNeighbourField()
                );

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = X_nei[facei] - X_[bCellID];

                    Ugrad[bCellID] +=
                        Ainv_[bCellID] & (U_nei[facei]-U[bCellID])*d;
                }
            }
        }

        Ugrad.correctBoundaryConditions();
    }

    tvf.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::gradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();

    UgradX = gradientSchemes::gradient(U.component(0));
    UgradY = gradientSchemes::gradient(U.component(1));
    UgradZ = gradientSchemes::gradient(U.component(2));

    tmp<GeometricField<tensor, fvPatchField, volMesh> > ttf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> Ugrad = ttf();

    forAll(mesh_.cells(), cellID)
    {
        Ugrad[cellID] = tensor(UgradX[cellID], UgradY[cellID], UgradZ[cellID]);
    }

    if( Pstream::parRun() )
    {
        Ugrad.correctBoundaryConditions();
    }

    tvf.clear();
    ttf.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
    const GeometricField<tensor, fvPatchField, volMesh>& U,
    GeometricField<tensor, fvPatchField, volMesh>& UgradX,
    GeometricField<tensor, fvPatchField, volMesh>& UgradY,
    GeometricField<tensor, fvPatchField, volMesh>& UgradZ
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    if (Pstream::parRun())
    {
        Ux.correctBoundaryConditions();
        Uy.correctBoundaryConditions();
        Uz.correctBoundaryConditions();
    }

    UgradX = gradientSchemes::gradient(Ux);
    UgradY = gradientSchemes::gradient(Uy);
    UgradZ = gradientSchemes::gradient(Uz);

    tvf.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::localGradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Unei
) const
{
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tvft
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> Ugrad = tvft();

    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];
        const vector& dOwn = XF_[faceID] - X_[ownID];
        const vector& dNei = XF_[faceID] - X_[neiID];

        UgradX[ownID] += AinvLocal_[ownID] & ((Unei[faceID].x()-U[ownID].x())*dOwn);
        UgradY[ownID] += AinvLocal_[ownID] & ((Unei[faceID].y()-U[ownID].y())*dOwn);
        UgradZ[ownID] += AinvLocal_[ownID] & ((Unei[faceID].z()-U[ownID].z())*dOwn);

        UgradX[neiID] += AinvLocal_[neiID] & ((Unei[faceID].x()-U[neiID].x())*dNei);
        UgradY[neiID] += AinvLocal_[neiID] & ((Unei[faceID].y()-U[neiID].y())*dNei);
        UgradZ[neiID] += AinvLocal_[neiID] & ((Unei[faceID].z()-U[neiID].z())*dNei);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = XF_.boundaryField()[patchID][facei] - X_[bCellID];

            UgradX[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].x() - U[bCellID].x())*d);

            UgradY[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].y() - U[bCellID].y())*d);

            UgradZ[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].z() - U[bCellID].z())*d);

            if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            {
                const label& faceID =
                    mesh_.boundary()[patchID].start() + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei];
                    vector d = XN_[nodeID] - X_[bCellID];

                    UgradX[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].x() - U[bCellID].x())*d);

                    UgradY[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].y() - U[bCellID].y())*d);

                    UgradZ[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].z() - U[bCellID].z())*d);

                    for (int i=0; i<7; i++)
                    {
                        d =
                            ((((i+1)*XN_[nodeID])
                          + ((7-i)*XF_.boundaryField()[patchID][facei]))/8.0)
                          - X_[bCellID];

                        UgradX[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].x() - U[bCellID].x())*d);

                        UgradY[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].y() - U[bCellID].y())*d);

                        UgradZ[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].z() - U[bCellID].z())*d);
                    }
                }
            }
        }
    }

    forAll(mesh_.cells(), cellID)
    {
        Ugrad[cellID] = tensor(UgradX[cellID], UgradY[cellID], UgradZ[cellID]);
    }

    tvf.clear();
    tvft.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Um,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Up
)
{
    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];

        Um[faceID] = U[ownID] + (Ugrad[ownID] & (XF_[faceID] - X_[ownID]));
        Up[faceID]  = U[neiID] + (Ugrad[neiID] & (XF_[faceID] - X_[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        forAll(mesh_.boundaryMesh()[patchID],facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            U.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + ( Ugrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            Um.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + ( Ugrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Um,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Up
)
{
    forAll(own_, faceID)
    {
        const label& ownID = own_[faceID];
        const label& neiID = nei_[faceID];

        Um[faceID] = U[ownID] + (Ugrad[ownID] & (XF_[faceID] - X_[ownID]));
        Up[faceID] = U[neiID] + (Ugrad[neiID] & (XF_[faceID] - X_[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            Um.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + (Ugrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<tensor, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Um,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Up
)
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            U.dimensions()
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "reconstruct("+Um.name()+')',
                Um.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Um.dimensions()
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh> UmX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmZ = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpZ = tsf();

    op.decomposeTensor(Um, UmX, UmY, UmZ);
    op.decomposeTensor(Up, UpX, UpY, UpZ);

    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ);

    forAll(own_, faceID)
    {
        Um[faceID] = tensor(UmX[faceID], UmY[faceID], UmZ[faceID]);
        Up[faceID] = tensor(UpX[faceID], UpY[faceID], UpZ[faceID]);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            const vector& reconsX =
                Ux[bCellID] + (UxGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            const vector& reconsY =
                Uy[bCellID] + (UyGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            const vector& reconsZ =
                Uz[bCellID] + (UzGrad[bCellID]
              & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));

            U.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);

            Um.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
        }
    }

    tvf.clear();
    tsf.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
