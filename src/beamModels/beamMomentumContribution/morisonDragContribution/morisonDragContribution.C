/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "morisonDragContribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(morisonDragContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, morisonDragContribution, dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


morisonDragContribution::morisonDragContribution
(
    const word& name,
    const dictionary& dict
)
:
    beamMomentumContribution
    (
        name,
        dict
    ),
    beamMomentumContribDict_(dict.subDict(name + "Coeffs")),
    Cdn_
    (
        readScalar(beamMomentumContribDict_.lookup("Cdn"))
    ),
    Cdt_
    (
        readScalar(beamMomentumContribDict_.lookup("Cdt"))
    ),
    rhoFluid_
    (
        readScalar(beamMomentumContribDict_.lookup("rhoFluid"))
    )
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<scalarSquareMatrix>> morisonDragContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<Field<scalarSquareMatrix>> tresult
    (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
    );
   // Field<scalarSquareMatrix>& result = tresult.ref();

    return tresult;
}


tmp<vectorField> morisonDragContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = bm.solutionW().mesh();

    // Flag to check the ddtSchemeName
    const word ddtSchemeName
    (
        mesh.ddtSchemes().found("ddt(W)")
        ?
        (mesh.ddtSchemes().lookup("ddt(W)"))
        :
        (mesh.ddtSchemes().lookup("default"))
    );

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    //- Drag forces due to Morison's Equation
    if (ddtSchemeName != "steadyState")
    {
        // Unit tangents to beam centreline at beam CV cell-centres
        const vectorField dRdScellHat(unitTangents(bm));

        // Tangential component of velocity vector
        vectorField Ut
        (
            (
                (U.internalField() & dRdScellHat)
                *dRdScellHat
            )
        );

        vectorField UtHat (Ut/(mag(Ut) + SMALL));

        // Normal component of velocity vector
        vectorField Un
        (
            (
                U.internalField()
                - (
                    (U.internalField() & dRdScellHat)
                    *dRdScellHat
                )
            )
        );

        vectorField UnHat (Un/(mag(Un) + SMALL));

        // Beam Radius and Length
        const scalar R = bm.R();
        const volScalarField& L = bm.L();

        // Scalar values of drag force (normal and tangential)
        const scalarField Fdn(rhoFluid_*Cdn_*R*L*(Un & Un));
        const scalarField Fdt(rhoFluid_*Cdt_*R*L*(Ut & Ut));

        forAll(result, cellI)
        {
            result[cellI] =
                Fdn[cellI]*UnHat[cellI]
              + Fdt[cellI]*UtHat[cellI];
        }
    }
    else
    {
        FatalErrorInFunction
        << "Momentum contribution due to Morison drag - type: "
        << typeName
        << " for time scheme " << ddtSchemeName
        << " is not defined!!\n"
        << "Set d2dt2Scheme and ddtScheme in system/fvSchemes as "
        << "'Euler' or 'Newmark'"
        << "to include drag force contributions"
        << abort(FatalError);
    }

    return tresult;
}


tmp<vectorField> morisonDragContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    // vectorField& result = tresult.ref();

    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
