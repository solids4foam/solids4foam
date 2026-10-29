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

#include "morisonAddedMassContribution.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(morisonAddedMassContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, morisonAddedMassContribution, dictionary
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

word morisonAddedMassContribution::d2dt2SchemeName(const fvMesh& mesh)
{
    const word schemeName
    (
        mesh.d2dt2Schemes().found("d2dt2(W)")
        ?
        word(mesh.d2dt2Schemes().lookup("d2dt2(W)"))
        :
        word(mesh.d2dt2Schemes().lookup("default"))
    );

    if (schemeName != "Euler" && schemeName != "Newmark")
    {
        FatalErrorInFunction
            << "Momentum contribution due to Morison added mass - type: "
            << typeName
            << " for time scheme " << schemeName
            << " is not defined!!\n"
            << "Set d2dt2Scheme and ddtScheme in system/fvSchemes as "
            << "'Euler' or 'Newmark' "
            << "to include added-mass force contributions"
            << abort(FatalError);
    }

    return schemeName;
}


tmp<scalarField> morisonAddedMassContribution::fluidMass
(
    const beamModel& bm
) const
{
    // pi*R^2 is the displaced area, and the added mass is the same in all
    // directions normal to the beam, only for a circular cross-section
    // (see requiresCircularSection)
    const scalar R = bm.R();

    return tmp<scalarField>
    (
        new scalarField(rhoFluid_*M_PI*R*R*bm.L().internalField())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

morisonAddedMassContribution::morisonAddedMassContribution
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
    CMn_
    (
        readScalar(beamMomentumContribDict_.lookup("CMn"))
    ),
    CMt_
    (
        readScalar(beamMomentumContribDict_.lookup("CMt"))
    ),
    rhoFluid_
    (
        readScalar(beamMomentumContribDict_.lookup("rhoFluid"))
    )
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<scalarSquareMatrix>> morisonAddedMassContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = bm.solutionW().mesh();

    const scalar deltaT = mesh.time().deltaTValue();

    // Derivative of the acceleration with respect to the displacement, as
    // used for the implicit beam inertia term
    scalar dAccldW = 1/sqr(deltaT);

    if (d2dt2SchemeName(mesh) == "Newmark")
    {
        const scalar betaN =
            mesh.d2dt2Schemes().lookupOrDefault<scalar>("newmarkBeta(W)", 0.25);

        dAccldW /= betaN;
    }

    // Prepare the result
    tmp<Field<scalarSquareMatrix>> tresult
    (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
    );
    Field<scalarSquareMatrix>& result = tresult.ref();

    const vectorField tHat(unitTangents(bm));

    // Displaced fluid mass of each beam CV
    const scalarField mf(fluidMass(bm));

    forAll(result, cellI)
    {
        // Added-mass tensor: CMn normal to the beam, CMt along it
        const tensor tt(tHat[cellI]*tHat[cellI]);
        const tensor addedMass(mf[cellI]*(CMn_*(I - tt) + CMt_*tt));

        // Implicit added-mass force coefficient (same sign as the beam
        // inertia coefficient QRhoCoeff)
        const tensor coeff(-dAccldW*addedMass);

        result[cellI](0,0) = coeff.xx();
        result[cellI](0,1) = coeff.xy();
        result[cellI](0,2) = coeff.xz();

        result[cellI](1,0) = coeff.yx();
        result[cellI](1,1) = coeff.yy();
        result[cellI](1,2) = coeff.yz();

        result[cellI](2,0) = coeff.zx();
        result[cellI](2,1) = coeff.zy();
        result[cellI](2,2) = coeff.zz();
    }

    return tresult;
}


tmp<vectorField> morisonAddedMassContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = bm.solutionW().mesh();

    // Beam acceleration, evaluated as in the beam inertia term
    const vectorField acclb
    (
        d2dt2SchemeName(mesh) == "Newmark"
      ? Accl.internalField()
      : fvc::ddt(U)().internalField()
    );

    const vectorField tHat(unitTangents(bm));

    // Tangential and normal components of the beam acceleration
    const vectorField acclbt((acclb & tHat)*tHat);
    const vectorField acclbn(acclb - acclbt);

    // Displaced fluid mass of each beam CV
    const scalarField mf(fluidMass(bm));

    // Added-mass force, with the same sign as the beam inertia force QRho
    return tmp<vectorField>
    (
        new vectorField(mf*(CMn_*acclbn + CMt_*acclbt))
    );
}


tmp<vectorField> morisonAddedMassContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    // No added rotational inertia
    return tmp<vectorField>
    (
        new vectorField(bm.solutionW().mesh().nCells(), vector::zero)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
