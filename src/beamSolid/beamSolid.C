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

#include "beamSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(beamSolid, 0);
addToRunTimeSelectionTable(solidModel, beamSolid, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void beamSolid::updateDisplacements()
{
    D() = beam_().solutionW();
    DD() = D() - D().oldTime();
    pointD() = mesh().lookupObject<pointVectorField>("pointW");
    pointDD() = pointD() - pointD().oldTime();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

beamSolid::beamSolid(Time& runTime, const word& region)
:
    solidModel(typeName, runTime, region),
    beam_(beamModel::New(runTime, region))
{
    // The beam model always uses its own Newton-Raphson solver, configured
    // in constant/beamProperties
    if (solidModelDict().found("solutionAlgorithm"))
    {
        WarningInFunction
            << "'solutionAlgorithm' in " << typeName << "Coeffs is ignored: "
            << "the beam is solved with the " << beam_().type()
            << " solver set in constant/beamProperties" << endl;
    }

    // Beam fields carry the physical state. The inherited
    // continuum fields are interface storage, not beam stress results.
    D().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    pointD().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
    sigma().writeOpt() = IOobject::NO_WRITE;

    beam_().updateTotalFields();
    updateDisplacements();
    D().oldTime() = D();
    pointD().oldTime() = pointD();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volVectorField& beamSolid::U()
{
    return const_cast<volVectorField&>
    (
        mesh().lookupObject<volVectorField>("U")
    );
}


const volVectorField& beamSolid::U() const
{
    return mesh().lookupObject<volVectorField>("U");
}


bool beamSolid::evolve()
{
    // The beam model returns an initial residual, not a success flag.
    // Non-convergence is reported by the beam model as a fatal error.
    beam_().evolve();
    return true;
}


void beamSolid::updateTotalFields()
{
    beam_().updateTotalFields();
    updateDisplacements();
}


void beamSolid::writeFields(const Time& runTime)
{
    beam_().writeFields();
}


tmp<vectorField> beamSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    FatalErrorInFunction
        << "Use beam force/moment boundary conditions for beamSolid. "
        << "Continuum traction boundary conditions are not supported."
        << abort(FatalError);
    return tmp<vectorField>(new vectorField(patch.size(), vector::zero));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
