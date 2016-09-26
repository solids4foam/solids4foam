/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fvcGradf.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvcSnGrad.H"
#include "gaussGrad.H"
//#include "ggiFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, 
        fvsPatchField, 
        surfaceMesh
    >
> fGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvsPatchField, surfaceMesh> > tGrad
    (
        new GeometricField<GradType, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "grad" + vf.name() + "f",
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
    Field<GradType>& gradI = tGrad().internalField();

    const vectorField& points = mesh.points();
    const faceList& faces = mesh.faces();

    surfaceVectorField n = mesh.Sf()/mesh.magSf();
    const vectorField& nI = n.internalField();

    const Field<Type>& pfI = pf.internalField();

    forAll(gradI, faceI)
    {
        const face& curFace = faces[faceI];

        vector Rf = curFace.centre(points);

        scalar mag = curFace.mag(points);

        const edgeList curFaceEdges = curFace.edges();

        gradI[faceI] = pTraits<GradType>::zero;
        scalar faceArea = 0;

        forAll(curFaceEdges, edgeI)
        {
            const edge& curEdge = curFaceEdges[edgeI];

            // Projected edge vector
            vector e = curEdge.vec(points);
            e -= nI[faceI]*(nI[faceI]&e);

            // Edge length vector
            vector Le = (e^nI[faceI]);
            Le *= curFace.edgeDirection(curEdge);

            // Edge-centre field value
            Type fe = 
                0.5
               *(
                   pfI[curEdge.start()] 
                 + pfI[curEdge.end()]
                );

            // Gradient
            gradI[faceI] += Le*fe;

            // Area
            vector Re = curEdge.centre(points) - Rf;
            Re -= nI[faceI]*(nI[faceI]&Re);
            faceArea += (Le&Re);
        }
        
        faceArea /= 2.0;

        gradI[faceI] /= mag; // faceArea; // mag
    }

    forAll(tGrad().boundaryField(), patchI)
    {
        Field<GradType>& patchGrad = tGrad().boundaryField()[patchI];

        const vectorField& patchN = n.boundaryField()[patchI];

        forAll(patchGrad, faceI)
        {
            label globalFaceID = mesh.boundaryMesh()[patchI].start() + faceI;

            const face& curFace = mesh.faces()[globalFaceID];

            vector Rf = curFace.centre(points);

            scalar mag = curFace.mag(points);

            const edgeList curFaceEdges = curFace.edges();

            patchGrad[faceI] = pTraits<GradType>::zero;
            scalar faceArea = 0;

            forAll(curFaceEdges, edgeI)
            {
                const edge& curEdge = curFaceEdges[edgeI];

                // Projected edge vector
                vector e = curEdge.vec(points);
                e -= patchN[faceI]*(patchN[faceI]&e);

                // Edge length vector
                vector Le = (e^patchN[faceI]);
                Le *= curFace.edgeDirection(curEdge);

                // Edge-centre field value
                Type fe = 
                    0.5
                   *(
                       pfI[curEdge.start()] 
                     + pfI[curEdge.end()]
                    );

                // Gradient
                patchGrad[faceI] += Le*fe;

                // Area
                vector Re = curEdge.centre(points) - Rf;
                Re -= patchN[faceI]*(patchN[faceI]&Re);
                faceArea += (Le&Re);
            }

            faceArea /= 2.0;
        
            patchGrad[faceI] /= mag; //faceArea; //mag
        }
    }

    // JN/PC: we need to fix this with AMI: for now throw an error
    forAll(mesh.boundary(), patchI)
    {
        //if (mesh.boundary()[patchI].type() == ggiFvPatch::typeName)
        if (mesh.boundary()[patchI].type() == "ggi")
        {
            FatalErrorIn("fvc::fGrad(...)")
                << "AMI functionality yet to be included in the of30 branch"
                << abort(FatalError);

            //const ggiFvPatch& ggiPatch =
            //    refCast<const ggiFvPatch>(mesh.boundary()[patchI]);

            //if (!ggiPatch.master())
            //{
            //    Field<GradType>& slaveGrad = 
            //        tGrad().boundaryField()[patchI];
            //    const Field<GradType>& masterGrad = 
            //        tGrad().boundaryField()[ggiPatch.shadowIndex()];

            //    slaveGrad = ggiPatch.interpolate(masterGrad);
            //}
        }
    }

//     const GeometricField<GradType, fvPatchField, volMesh>& gradU = 
//         mesh.lookupObject<GeometricField<GradType, fvPatchField, volMesh> >
//         (
//             "grad(" + vf.name() + ")"
//         );
//     tGrad() = ((I - n*n)&fvc::interpolate(gradU));


    // Add normal component of the gradient
    tGrad() += n*fvc::snGrad(vf);

    return tGrad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, 
        fvsPatchField, 
        surfaceMesh
    >
> fsGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvsPatchField, surfaceMesh> > tGrad
    (
        new GeometricField<GradType, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "grad" + vf.name() + "f",
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
    Field<GradType>& gradI = tGrad().internalField();

    const vectorField& points = mesh.points();
    const faceList& faces = mesh.faces();

    surfaceVectorField n = mesh.Sf()/mesh.magSf();
    const vectorField& nI = n.internalField();

    const Field<Type>& pfI = pf.internalField();

    forAll(gradI, faceI)
    {
        const face& curFace = faces[faceI];

        vector Rf = curFace.centre(points);

        scalar mag = curFace.mag(points);

        const edgeList curFaceEdges = curFace.edges();

        gradI[faceI] = pTraits<GradType>::zero;
        scalar faceArea = 0;

        forAll(curFaceEdges, edgeI)
        {
            const edge& curEdge = curFaceEdges[edgeI];

            // Projected edge vector
            vector e = curEdge.vec(points);
            e -= nI[faceI]*(nI[faceI]&e);

            // Edge length vector
            vector Le = (e^nI[faceI]);
            Le *= curFace.edgeDirection(curEdge);

            // Edge-centre field value
            Type fe = 
                0.5
               *(
                   pfI[curEdge.start()] 
                 + pfI[curEdge.end()]
                );

            // Gradient
            gradI[faceI] += Le*fe;

            // Area
            vector Re = curEdge.centre(points) - Rf;
            Re -= nI[faceI]*(nI[faceI]&Re);
            faceArea += (Le&Re);
        }
        
        faceArea /= 2.0;

        gradI[faceI] /= mag; // faceArea; // mag
    }

    forAll(tGrad().boundaryField(), patchI)
    {
        Field<GradType>& patchGrad = tGrad().boundaryField()[patchI];

        const vectorField& patchN = n.boundaryField()[patchI];

        forAll(patchGrad, faceI)
        {
            label globalFaceID = mesh.boundaryMesh()[patchI].start() + faceI;

            const face& curFace = mesh.faces()[globalFaceID];

            vector Rf = curFace.centre(points);

            scalar mag = curFace.mag(points);

            const edgeList curFaceEdges = curFace.edges();

            patchGrad[faceI] = pTraits<GradType>::zero;
            scalar faceArea = 0;

            forAll(curFaceEdges, edgeI)
            {
                const edge& curEdge = curFaceEdges[edgeI];

                // Projected edge vector
                vector e = curEdge.vec(points);
                e -= patchN[faceI]*(patchN[faceI]&e);

                // Edge length vector
                vector Le = (e^patchN[faceI]);
                Le *= curFace.edgeDirection(curEdge);

                // Edge-centre field value
                Type fe = 
                    0.5
                   *(
                       pfI[curEdge.start()] 
                     + pfI[curEdge.end()]
                    );

                // Gradient
                patchGrad[faceI] += Le*fe;

                // Area
                vector Re = curEdge.centre(points) - Rf;
                Re -= patchN[faceI]*(patchN[faceI]&Re);
                faceArea += (Le&Re);
            }

            faceArea /= 2.0;
        
            patchGrad[faceI] /= mag; //faceArea; //mag
        }
    }

    // JN/PC: we need to fix this with AMI: for now throw an error
    forAll(mesh.boundary(), patchI)
    {
        //if (mesh.boundary()[patchI].type() == ggiFvPatch::typeName)
        if (mesh.boundary()[patchI].type() == "ggi")
        {
            FatalErrorIn("fvc::fGrad(...)")
                << "AMI functionality yet to be included in the of30 branch"
                << abort(FatalError);

            //const ggiFvPatch& ggiPatch =
            //    refCast<const ggiFvPatch>(mesh.boundary()[patchI]);

            //if (!ggiPatch.master())
            //{
            //    Field<GradType>& slaveGrad = 
            //        tGrad().boundaryField()[patchI];
            //    const Field<GradType>& masterGrad = 
            //        tGrad().boundaryField()[ggiPatch.shadowIndex()];

            //    slaveGrad = ggiPatch.interpolate(masterGrad);
            //}
        }
    }

    return tGrad;
}


template<class Type, template<class> class FaceList>
tmp<Field<typename outerProduct<vector, Type>::type> > fGrad
(
    const PrimitivePatch<face, FaceList, const pointField&>& patch,
    const Field<Type>& ppf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<Field<GradType> > tGrad
    (
        new Field<GradType>
        (
            patch.size(),
            pTraits<GradType>::zero
        )
    );
    Field<GradType>& grad = tGrad();

    const vectorField& points = patch.localPoints();
    const faceList& faces = patch.localFaces();

    forAll(grad, faceI)
    {
        const face& curFace = faces[faceI];

        vector n = curFace.normal(points);
        n /= mag(n);

        vector Rf = curFace.centre(points);

        scalar mag = curFace.mag(points);

        const edgeList curFaceEdges = curFace.edges();

        scalar faceArea = 0;

        forAll(curFaceEdges, edgeI)
        {
            const edge& curEdge = curFaceEdges[edgeI];

            // Projected edge vector
            vector e = curEdge.vec(points);
            e -= n*(n&e);

            // Edge length vector
            vector Le = (e^n);
            Le *= curFace.edgeDirection(curEdge);

            // Edge-centre displacement
            Type fe = 
                0.5
               *(
                   ppf[curEdge.start()]
                 + ppf[curEdge.end()]
                );

            // Gradient
            grad[faceI] += Le*fe;

            // Area
            vector Re = curEdge.centre(points) - Rf;
            Re -= n*(n&Re);
            faceArea += (Le&Re);
        }

        faceArea /= 2.0;
        
        grad[faceI] /= mag; //faceArea;
    }

    return tGrad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, 
        fvPatchField, 
        volMesh
    >
> grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "grad(" + vf.name() + ")",
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
            Field<Type> ppf = pf.boundaryField()[patchI].patchInternalField();

            tGrad().boundaryField()[patchI] = 
                fGrad(mesh.boundaryMesh()[patchI], ppf);
        }
    }

//     tGrad().correctBoundaryConditions();

    // Normal gradient
    fv::gaussGrad<Type>(mesh).correctBoundaryConditions(vf, tGrad());

    return tGrad;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
