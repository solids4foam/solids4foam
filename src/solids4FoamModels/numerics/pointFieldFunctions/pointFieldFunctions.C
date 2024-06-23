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

#include "pointFieldFunctions.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<pointSymmTensorField> symm(const pointTensorField& ptr)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointSymmTensorField> tresult
    (
        new pointSymmTensorField
        (
            IOobject
            (
                "symm(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("0", ptr.dimensions(), symmTensor::zero)
        )
    );

    // Set the result field to be symm(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = symm(ptr.primitiveField());
    tresult.ref().correctBoundaryConditions();
#else
    tresult().internalField() = symm(ptr.internalField());
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


tmp<pointSymmTensorField> dev(const pointSymmTensorField& ptr)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointSymmTensorField> tresult
    (
        new pointSymmTensorField
        (
            IOobject
            (
                "dev(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("0", ptr.dimensions(), symmTensor::zero)
        )
    );

    // Set the result field to be symm(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = dev(ptr.primitiveField());
#else
    tresult().internalField() = dev(ptr.internalField());
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


tmp<pointScalarField> tr(const pointSymmTensorField& ptr)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointScalarField> tresult
    (
        new pointScalarField
        (
            IOobject
            (
                "tr(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", ptr.dimensions(), 0.0)
        )
    );

    // Set the result field to be symm(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = tr(ptr.primitiveField());
    tresult.ref().correctBoundaryConditions();
#else
    tresult().internalField() = tr(ptr.internalField());
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


tmp<pointTensorField> cof(const pointTensorField& ptr)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointTensorField> tresult
    (
        new pointTensorField
        (
            IOobject
            (
                "cof(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedTensor
            (
                "zero", pow(ptr.dimensions(), tensor::dim - 1), tensor::zero
            )
        )
    );

    // Set the result field to be cof(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = cof(ptr.primitiveField());
    tresult.ref().correctBoundaryConditions();
#else
    tresult().internalField() = cof(ptr.internalField());
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


tmp<pointScalarField> sqr(const pointScalarField& ptr)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointScalarField> tresult
    (
        new pointScalarField
        (
            IOobject
            (
                "sqr(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", pow(ptr.dimensions(), 2), 0.0)
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = sqr(ptr.primitiveField());
    tresult.ref().correctBoundaryConditions();
#else
    tresult().internalField() = sqr(ptr.internalField());
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


tmp<pointScalarField> pow(const pointScalarField& ptr, const scalar& exponent)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointScalarField> tresult
    (
        new pointScalarField
        (
            IOobject
            (
                "pow(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", pow(ptr.dimensions(), exponent), 0.0)
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = pow(ptr.primitiveField(), exponent);
    tresult.ref().correctBoundaryConditions();
#else
    tresult().internalField() = pow(ptr.internalField(), exponent);
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


tmp<pointScalarField> sqrt(const pointScalarField& ptr)
{
    const pointMesh& pMesh =  ptr.mesh();

    // Prepare the temporary field
    tmp<pointScalarField> tresult
    (
        new pointScalarField
        (
            IOobject
            (
                "sqrt(" + ptr.name() + ")",
                pMesh.mesh().time().timeName(),
                pMesh.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", sqrt(ptr.dimensions()), 0.0)
        )
    );

    // Set the result field to be sqrt(ptr)
#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().primitiveFieldRef() = sqrt(ptr.primitiveField());
    tresult.ref().correctBoundaryConditions();
#else
    tresult().internalField() = sqrt(ptr.internalField());
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}

}

// ************************************************************************* //
