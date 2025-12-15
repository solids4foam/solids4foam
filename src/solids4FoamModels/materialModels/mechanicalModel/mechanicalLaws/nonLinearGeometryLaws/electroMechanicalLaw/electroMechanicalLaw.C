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

#include "electroMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, electroMechanicalLaw, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::electroMechanicalLaw::electroMechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    passiveMechLawPtr_
    (
        mechanicalLaw::NewNonLinGeomMechLaw
        (
            word(dict.subDict("passiveMechanicalLaw").lookup("type")),
            mesh,
            dict.subDict("passiveMechanicalLaw"),
            nonLinGeom
        )
    ),
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime")))
{
    if (rampTime_ < 0.0)
    {
        FatalErrorIn("electroMechanicalLaw::electroMechanicalLaw(...)")
            << "rampTime should be greater than or equal to zero"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroMechanicalLaw::~electroMechanicalLaw()
{}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::impK() const
{
    return passiveMechLawPtr_->impK();
}


#ifdef OPENFOAM_NOT_EXTEND
Foam::tmp<Foam::Field<Foam::scalarSquareMatrix>>
Foam::electroMechanicalLaw::materialTangentField() const
{
    return passiveMechLawPtr_->materialTangentField();
}
#endif


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::bulkModulus() const
{
    return passiveMechLawPtr_->bulkModulus();
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::shearModulus() const
{
    return passiveMechLawPtr_->shearModulus();
}


void Foam::electroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Lookup the fibre directions
    // How best should we do this to avoid duplicating the fibre field?
    // For now, let's hard-code in the field name
    const volSymmTensorField f0f0 =
        mesh().lookupObject<volSymmTensorField>("f0f0");

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // For now, we will assume a constant active stress
    // The next step will be to include an active-stress model to convert
    // muscle activation to fibre tension

    // Calculate current value of Ta
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    // Add active stress to the passive stress
    // Note that the active stress is converted from a 2nd Piola-Kirchhoff
    // stress to a Cauchy stress
    // sigma += J*symm(F & (currentTa*f0f0) & F.T());
    sigma += symm(F & (currentTa*f0f0) & F.T())/J;
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Lookup the fibre directions
    // How best should we do this to avoid duplicating the fibre field?
    // For now, let's hard-code in the field name
    const surfaceSymmTensorField f0f0 =
        mesh().lookupObject<surfaceSymmTensorField>("f0f0f");

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const surfaceTensorField& F = this->Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    // For now, we will assume a constant active stress
    // The next step will be to include an active-stress model to convert
    // muscle activation to fibre tension

    // Calculate current value of Ta
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    // Add active stress to the passive stress
    // Note that the active stress is converted from a 2nd Piola-Kirchhoff
    // stress to a Cauchy stress
    sigma += J*symm(F & (currentTa*f0f0) & F.T());
}

void Foam::electroMechanicalLaw::correct
(
    List<List<symmTensor>>& sigmaQuad,
    const List<List<tensor>>& gradDQuad
)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigmaQuad, gradDQuad);

    // Calculate current value of Ta
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    // Lookup the fibre directions
    // How best should we do this to avoid duplicating the fibre field?
    // For now, let's hard-code in the field name
    //    const surfaceSymmTensorField f0f0f =
    //    mesh().lookupObject<surfaceSymmTensorField>("f0f0f");
    symmTensor f0f0f(0, 0, 0, 0, 0, 1.0);

    // Initialise outside loop for efficiency
    tensor F = tensor::zero;
    scalar J = 0.0;

    forAll(sigmaQuad, faceI)
    {
        List<symmTensor>& faceSigmaQuad = sigmaQuad[faceI];
        const List<tensor>& faceGradDQuad = gradDQuad[faceI];

        forAll(faceSigmaQuad, qpI)
        {
            F = I + faceGradDQuad[qpI].T();
            J = det(F);
            // if (J < SMALL)
            // {
	    // 	WarningInFunction
            //         << "Unphysical Jacobian J = " << J << " at face " << faceI
            //         << ", qp " << qpI << "." << endl;
	    // }
	    sigmaQuad[faceI][qpI] = J * symm(F & (currentTa.value()*f0f0f) & F.T());
        }
    }
}


// ************************************************************************* //
