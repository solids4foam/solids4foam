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

#include "neoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::neoHookeanElastic::calculateStress
(
 surfaceSymmTensorField& sigma,
 const surfaceTensorField& gradD
)
{
	//Calculate F
	Ff() = I + gradD.T();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField bEbar(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(mu_*dev(bEbar));

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElastic::neoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0)
{
    // Read mechanical properties
    if
    (
        dict.found("E") && dict.found("nu")
     && !dict.found("mu") && !dict.found("K")
    )
    {
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        mu_ = (E/(2.0*(1.0 + nu)));

        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if
    (
        dict.found("mu") && dict.found("K")
     && !dict.found("E") && !dict.found("nu")
    )
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn(type())
            << "Either E and nu or mu and K should be specified"
            << abort(FatalError);
    }

    // Store old F
    F().storeOldTime();
    Ff().storeOldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElastic::~neoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
     new volScalarField
     (
      IOobject
      (
       "impK",
       mesh().time().timeName(),
       mesh(),
       IOobject::NO_READ,
       IOobject::NO_WRITE
       ),
      mesh(),
      K_
      )
     );
}

Foam::tmp<Foam::Field<Foam::RectangularMatrix<Foam::scalar>>>
Foam::neoHookeanElastic::materialTangentField() const
{
    // Prepare tmp field
    tmp<Field<Foam::RectangularMatrix<Foam::scalar>>> tresult
    (
        new Field<Foam::RectangularMatrix<Foam::scalar>>(mesh().nFaces(), Foam::RectangularMatrix<scalar>(6,9,0))
    );
#ifdef OPENFOAM_NOT_EXTEND
    Field<Foam::RectangularMatrix<Foam::scalar>>& result = tresult.ref();
#else
    Field<Foam::RectangularMatrix<Foam::scalar>>& result = tresult();
#endif

    // Calculate tangent field
    const Switch numericalTangent(dict().lookup("numericalTangent"));
    if (numericalTangent)
    {
        // Lookup current stress and store it as the reference
        const surfaceSymmTensorField& sigmaRef =
            mesh().lookupObject<surfaceSymmTensorField>("sigmaf");
        // Lookup gradient of displacement
        const surfaceTensorField& gradDRef =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Create fields to be used for perturbations
        surfaceSymmTensorField sigmaPerturb("sigmaPerturb", sigmaRef);
        surfaceTensorField gradDPerturb("gradDPerturb", gradDRef);

        // Small number used for perturbations
        const scalar eps(readScalar(dict().lookup("tangentEps")));

        // For each component of gradD, sequentially apply a perturbation and
        // then calculate the resulting sigma
        for (label cmptI = 0; cmptI < tensor::nComponents; cmptI++)
        {
			// Reset gradDPerturb and multiply by 1.0 to avoid it being removed
			// from the object registry
			gradDPerturb = 1.0*gradDRef;

            // Perturb this component of gradD and calculate FPerturb
            gradDPerturb.replace(cmptI, gradDRef.component(cmptI) + eps);

            // Calculate perturbed stress
            const_cast<neoHookeanElastic&>(*this).calculateStress(sigmaPerturb, gradDPerturb);

            // Calculate tangent component
            const surfaceSymmTensorField tangCmpt((sigmaPerturb - sigmaRef)/eps);
            const symmTensorField& tangCmptI = tangCmpt.internalField();

            // Insert tangent component
            forAll(tangCmptI, faceI)
            {
                if (cmptI == tensor::XX)
                {
                    result[faceI](0,0) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,0) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,0) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,0) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,0) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,0) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::XY)
                {
                    result[faceI](0,1) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,1) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,1) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,1) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,1) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,1) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::XZ)
                {
                    result[faceI](0,2) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,2) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,2) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,2) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,2) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,2) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::YX)
                {
                    result[faceI](0,3) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,3) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,3) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,3) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,3) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,3) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::YY)
                {
                    result[faceI](0,4) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,4) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,4) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,4) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,4) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,4) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::YZ)
                {
                    result[faceI](0,5) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,5) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,5) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,5) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,5) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,5) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::ZX)
                {
                    result[faceI](0,6) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,6) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,6) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,6) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,6) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,6) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == tensor::ZY)
                {
                    result[faceI](0,7) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,7) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,7) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,7) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,7) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,7) = tangCmptI[faceI][symmTensor::XZ];
                }
                else // if (cmptI == tensor::ZZ)
                {
                    result[faceI](0,8) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,8) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,8) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,8) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,8) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,8) = tangCmptI[faceI][symmTensor::XZ];
                }
            }

            forAll(tangCmpt.boundaryField(), patchI)
            {
                const symmTensorField& tangCmptP =
                    tangCmpt.boundaryField()[patchI];
                const label start = mesh().boundaryMesh()[patchI].start();

                forAll(tangCmptP, fI)
                {
                    const label faceID = start + fI;

		            if (cmptI == tensor::XX)
		            {
		                result[faceID](0,0) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,0) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,0) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,0) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,0) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,0) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::XY)
		            {
		                result[faceID](0,1) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,1) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,1) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,1) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,1) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,1) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::XZ)
		            {
		                result[faceID](0,2) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,2) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,2) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,2) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,2) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,2) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::YX)
		            {
		                result[faceID](0,3) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,3) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,3) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,3) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,3) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,3) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::YY)
		            {
		                result[faceID](0,4) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,4) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,4) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,4) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,4) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,4) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::YZ)
		            {
		                result[faceID](0,5) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,5) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,5) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,5) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,5) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,5) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::ZX)
		            {
		                result[faceID](0,6) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,6) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,6) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,6) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,6) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,6) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else if (cmptI == tensor::ZY)
		            {
		                result[faceID](0,7) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,7) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,7) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,7) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,7) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,7) = tangCmptI[fI][symmTensor::XZ];
		            }
		            else // if (cmptI == tensor::ZZ)
		            {
		                result[faceID](0,8) = tangCmptI[fI][symmTensor::XX];
		                result[faceID](1,8) = tangCmptI[fI][symmTensor::YY];
		                result[faceID](2,8) = tangCmptI[fI][symmTensor::ZZ];
		                result[faceID](3,8) = tangCmptI[fI][symmTensor::XY];
		                result[faceID](4,8) = tangCmptI[fI][symmTensor::YZ];
		                result[faceID](5,8) = tangCmptI[fI][symmTensor::XZ];
		            }
		        }
            }
        }
    }
    else // Analytical tangent
    {

        notImplemented("Analytical tangent not implemented");

    }

    return tresult;
}


void Foam::neoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField bEbar(pow(J, -2.0/3.0)*symm(F() & F().T()));

    // Calculate the deviatoric stress
    const volSymmTensorField s(mu_*dev(bEbar));

    // Update the hydrostatic stress
    updateSigmaHyd
    (
        0.5*K()*(pow(J, 2.0) - 1.0),
        (4.0/3.0)*mu_ + K_
    );

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(sigmaHyd()*I + s);
}


void Foam::neoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField bEbar(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(mu_*dev(bEbar));

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


void Foam::neoHookeanElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
