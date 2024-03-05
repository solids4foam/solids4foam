/*---------------------------------------------------------------------------* \
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

#include "pointGaussLeastSquaresGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "processorFvPatch.H"
#include "ggiFvPatch.H"
#include "regionCoupleFvPatch.H"
#include "fvc.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
pointGaussLeastSquaresGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                vf.dimensions()/dimLength,
                pTraits<GradType>::zero
            )
        )
    );

    Field<GradType>& iGrad = tGrad().internalField();

    // Lookup point field from the solver

    const GeometricField<Type, pointPatchField, pointMesh>& pf =
        mesh.lookupObject< GeometricField<Type, pointPatchField, pointMesh> >
        (
            "point" + vf.name()
        );

    const vectorField& points = mesh.points();

    const faceList& faces = mesh.faces();

    const Field<Type>& pfI = pf.internalField();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    scalarField V(iGrad.size(), 0.0);

    forAll(owner, faceI)
    {
        const face& curFace = faces[faceI];

        // If the face is a triangle, do a direct calculation
        if (curFace.size() == 3)
        {
            GradType SF =
                curFace.normal(points)*curFace.average(points, pfI);

            iGrad[owner[faceI]] += SF;
            iGrad[neighbour[faceI]] -= SF;

            scalar SR = (curFace.normal(points)&curFace.centre(points));

            V[owner[faceI]] += SR;
            V[neighbour[faceI]] -= SR;
        }
        else
        {
            label nPoints = curFace.size();

            point centrePoint = point::zero;
            Type cf = pTraits<Type>::zero;

            for (register label pI=0; pI<nPoints; pI++)
            {
                centrePoint += points[curFace[pI]];
                cf += pfI[curFace[pI]];
            }

            centrePoint /= nPoints;
            cf /= nPoints;

            for (register label pI=0; pI<nPoints; pI++)
            {
                // Calculate triangle centre field value
                Type ttcf =
                (
                    pfI[curFace[pI]]
                  + pfI[curFace[(pI + 1) % nPoints]]
                  + cf
                );
                ttcf /= 3.0;

                // Calculate triangle area
                vector St =
                    (
                        (points[curFace[pI]] - centrePoint)
                      ^ (
                            points[curFace[(pI + 1) % nPoints]]
                          - centrePoint
                        )
                    );
                St /= 2.0;

                // Calculate triangle centre
                vector Ct =
                    (
                        centrePoint
                      + points[curFace[pI]]
                      + points[curFace[(pI + 1) % nPoints]]
                    )/3;


                iGrad[owner[faceI]] += St*ttcf;
                iGrad[neighbour[faceI]] -= St*ttcf;

                V[owner[faceI]] += (St&Ct);
                V[neighbour[faceI]] -= (St&Ct);
            }
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundaryMesh()[patchI].faceCells();

        forAll(mesh.boundaryMesh()[patchI], faceI)
        {
            label globalFaceID =
                mesh.boundaryMesh()[patchI].start() + faceI;

            const face& curFace = faces[globalFaceID];

            // If the face is a triangle, do a direct calculation
            if (curFace.size() == 3)
            {
                iGrad[pFaceCells[faceI]] +=
                    curFace.normal(points)*curFace.average(points, pfI);

                V[pFaceCells[faceI]] +=
                    (curFace.normal(points)&curFace.centre(points));
            }
            else
            {
                label nPoints = curFace.size();

                point centrePoint = point::zero;
                Type cf = pTraits<Type>::zero;

                for (register label pI=0; pI<nPoints; pI++)
                {
                    centrePoint += points[curFace[pI]];
                    cf += pfI[curFace[pI]];
                }

                centrePoint /= nPoints;
                cf /= nPoints;

                for (register label pI=0; pI<nPoints; pI++)
                {
                    // Calculate triangle centre field value
                    Type ttcf =
                    (
                        pfI[curFace[pI]]
                      + pfI[curFace[(pI + 1) % nPoints]]
                      + cf
                    );
                    ttcf /= 3.0;

                    // Calculate triangle area
                    vector St =
                        (
                            (points[curFace[pI]] - centrePoint)
                          ^ (
                                points[curFace[(pI + 1) % nPoints]]
                              - centrePoint
                            )
                        );
                    St /= 2.0;

                    // Calculate triangle centre
                    vector Ct =
                        (
                            centrePoint
                          + points[curFace[pI]]
                          + points[curFace[(pI + 1) % nPoints]]
                        )/3;

                    iGrad[pFaceCells[faceI]] += St*ttcf;

                    V[pFaceCells[faceI]] += (St&Ct);
                }
            }
        }
    }

    V /= 3;

    iGrad /= V;
//     iGrad /= mesh.V();

//     iGrad = fv::gaussGrad<vector>(mesh).grad(vf)().internalField();

    // Calculate boundary gradient
    forAll(mesh.boundary(), patchI)
    {
        if
        (
            mesh.boundary()[patchI].size()
         && !vf.boundaryField()[patchI].coupled()
        )
        {
            // Info << mesh.boundary()[patchI].name() << ", "
            //      << mesh.boundary()[patchI].type() << endl;

            Field<Type> ppf = pf.boundaryField()[patchI].patchInternalField();

            tGrad().boundaryField()[patchI] =
                fvc::fGrad(mesh.boundaryMesh()[patchI], ppf);
        }
    }

    //     tGrad().correctBoundaryConditions();

    // Normal gradient
    fv::gaussGrad<Type>(mesh).correctBoundaryConditions(vf, tGrad());

    // Correct gradient at ggi patches
    forAll(mesh.boundary(), patchI)
    {
        if
        (
            isA<ggiFvPatch>(mesh.boundary()[patchI])
         || isA<regionCoupleFvPatch>(mesh.boundary()[patchI])
        )
        {
            Field<Type> ppf = pf.boundaryField()[patchI].patchInternalField();

            tGrad().boundaryField()[patchI] =
                fvc::fGrad(mesh.boundaryMesh()[patchI], ppf)
              + mesh.boundary()[patchI].nf()
               *vf.boundaryField()[patchI].snGrad();
        }
    }

    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].type() == ggiFvPatch::typeName)
        {
            const ggiFvPatch& ggiPatch =
                refCast<const ggiFvPatch>(mesh.boundary()[patchI]);

            if (!ggiPatch.master())
            {
                Field<GradType>& slaveGrad =
                    tGrad().boundaryField()[patchI];

                const Field<GradType>& masterGrad =
                    tGrad().boundaryField()[ggiPatch.shadowIndex()];

                slaveGrad = ggiPatch.interpolate(masterGrad);
            }
        }
    }

    return tGrad;
}


// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector, Type>::type, fvsPatchField, surfaceMesh
//     >
// >
// pointGaussLeastSquaresGrad<Type>::fGrad
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vsf
// ) const
// {
//     typedef typename outerProduct<vector, Type>::type GradType;

//     const fvMesh& mesh = vsf.mesh();

//     tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
//     (
//         new GeometricField<GradType, fvPatchField, volMesh>
//         (
//             IOobject
//             (
//                 "grad("+vsf.name()+')',
//                 vsf.instance(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh,
//             dimensioned<GradType>
//             (
//                 "zero",
//                 vsf.dimensions()/dimLength,
//                 pTraits<GradType>::zero
//             ),
//             zeroGradientFvPatchField<GradType>::typeName
//         )
//     );
//     //GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad();


//     return tlsGrad;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
