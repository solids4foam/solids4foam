/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvcDivGaussPoints.H"
#include "fvMesh.H"
#include "fvcSurfaceIntegrateGaussQuadIntegrate.H"

// We can remove schemes later and leave the one we use
#include "divScheme.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
divGaussPoints
(
    const List<List<symmTensor>>& sigmaGP,
    const List<List<scalar>>>& gpW
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            "divGaussPoints("+ssf.name()+')',
            fvc::surfaceGaussQuadIntegrate(sigamGP, gpW)
        )
    );
}


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
// )
// {
//     return tmp<GeometricField<Type, fvPatchField, volMesh>>
//     (
//         new GeometricField<Type, fvPatchField, volMesh>
//         (
//             "divGaussPoints("+ssf.name()+')',
//             fvc::surfaceIntegrate(ssf)
//         )
//     );
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints(fvc::divGaussPoints(tssf()));
//     tssf.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename innerProduct<vector, Type>::type, fvPatchField, volMesh
//     >
// >
// divGaussPoints
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vf,
//     const word& name
// )
// {
//     return fv::divGaussPointsScheme<Type>::New
//     (
//         vf.mesh(), vf.mesh().divScheme(name)
//     ).ref().fvcDivGaussPoints(vf);
// }


// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename innerProduct<vector, Type>::type, fvPatchField, volMesh
//     >
// >
// divGaussPoints
// (
//     const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvvf,
//     const word& name
// )
// {
//     typedef typename innerProduct<vector, Type>::type DivGaussPointsType;
//     tmp<GeometricField<DivGaussPointsType, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(tvvf(), name)
//     );
//     tvvf.clear();
//     return DivGaussPoints;
// }

// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename innerProduct<vector, Type>::type, fvPatchField, volMesh
//     >
// >
// divGaussPoints
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// )
// {
//     return fvc::divGaussPoints(vf, "divGaussPoints("+vf.name()+')');
// }


// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename innerProduct<vector, Type>::type, fvPatchField, volMesh
//     >
// >
// divGaussPoints
// (
//     const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvvf
// )
// {
//     typedef typename innerProduct<vector, Type>::type DivGaussPointsType;
//     tmp<GeometricField<DivGaussPointsType, fvPatchField, volMesh>> DivGaussPoints(fvc::divGaussPoints(tvvf()));
//     tvvf.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const surfaceScalarField& flux,
//     const GeometricField<Type, fvPatchField, volMesh>& vf,
//     const word& name
// )
// {
//     return fv::convectionScheme<Type>::New
//     (
//         vf.mesh(),
//         flux,
//         vf.mesh().divGaussPointsScheme(name)
//     ).ref().fvcDivGaussPoints(flux, vf);
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const tmp<surfaceScalarField>& tflux,
//     const GeometricField<Type, fvPatchField, volMesh>& vf,
//     const word& name
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(tflux(), vf, name)
//     );
//     tflux.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const surfaceScalarField& flux,
//     const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
//     const word& name
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(flux, tvf(), name)
//     );
//     tvf.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const tmp<surfaceScalarField>& tflux,
//     const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
//     const word& name
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(tflux(), tvf(), name)
//     );
//     tflux.clear();
//     tvf.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const surfaceScalarField& flux,
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// )
// {
//     return fvc::divGaussPoints
//     (
//         flux, vf, "divGaussPoints("+flux.name()+','+vf.name()+')'
//     );
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const tmp<surfaceScalarField>& tflux,
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(tflux(), vf)
//     );
//     tflux.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const surfaceScalarField& flux,
//     const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(flux, tvf())
//     );
//     tvf.clear();
//     return DivGaussPoints;
// }


// template<class Type>
// tmp<GeometricField<Type, fvPatchField, volMesh>>
// divGaussPoints
// (
//     const tmp<surfaceScalarField>& tflux,
//     const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
// )
// {
//     tmp<GeometricField<Type, fvPatchField, volMesh>> DivGaussPoints
//     (
//         fvc::divGaussPoints(tflux(), tvf())
//     );
//     tflux.clear();
//     tvf.clear();
//     return DivGaussPoints;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
