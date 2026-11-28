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

#include "logStrainMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "eig3.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(logStrainMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        logStrainMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::logStrainMechanicalConstitutiveLaw::trialElasticStrain
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawState& state,
    Field<symmTensor>& epsilonETrial,
    Field<tensor>& eigenVectors
) const
{
    const UIndirectList<tensor>& F = kin.F();
    const UIndirectList<tensor>& Finv0 = kin.Finv0();
    const UIndirectList<scalar>& J = kin.J();

    const Field<symmTensor>& bE0 = state.getSymmTensorField0("bE");

    forAll(F, i)
    {
        // An inverted point has no logarithmic strain, and f & bE0 & f^T
        // would hide it, being positive definite for any invertible f
        if (J[i] <= 0)
        {
            FatalErrorInFunction
                << "Non-positive Jacobian J = " << J[i]
                << " at integration point " << i << ": the deformation has "
                << "inverted the material, which has no logarithmic strain"
                << exit(FatalError);
        }

        // Relative deformation gradient over the step
        const tensor relF(F[i] & Finv0[i]);

        // Elastic left Cauchy-Green tensor at the start of the step, pushed
        // forward by the step's deformation as if it were elastic
        const symmTensor bETrial(symm(relF & bE0[i] & relF.T()));

        // Its logarithm, keeping the eigenvectors for the update of bE.
        // eig3 rather than OpenFOAM's eigenVectors, which are unreliable for
        // repeated eigenvalues, as at the undeformed state. Each eigenvector
        // is a row of V
        tensor& V = eigenVectors[i];
        vector d(vector::zero);
        eig3::eigen_decomposition(bETrial, V, d);

        if (min(d.x(), min(d.y(), d.z())) <= 0)
        {
            FatalErrorInFunction
                << "The trial elastic left Cauchy-Green tensor " << bETrial
                << " at integration point " << i << " is not positive "
                << "definite: its eigenvalues are " << d
                << exit(FatalError);
        }

        epsilonETrial[i] =
            0.5
           *(
                Foam::log(d.x())*sqr(V.x())
              + Foam::log(d.y())*sqr(V.y())
              + Foam::log(d.z())*sqr(V.z())
            );
    }
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

Foam::symmTensor Foam::logStrainMechanicalConstitutiveLaw::logm
(
    const symmTensor& A
)
{
    // eig3 rather than OpenFOAM's eigenVectors, which are unreliable for
    // repeated eigenvalues, as at the undeformed state. Each eigenvector is
    // a row of V
    tensor V(tensor::zero);
    vector d(vector::zero);
    eig3::eigen_decomposition(A, V, d);

    if (min(d.x(), min(d.y(), d.z())) <= 0)
    {
        FatalErrorInFunction
            << "The logarithm of " << A << " is not defined: its "
            << "eigenvalues " << d << " are not all positive"
            << exit(FatalError);
    }

    return
        Foam::log(d.x())*sqr(V.x())
      + Foam::log(d.y())*sqr(V.y())
      + Foam::log(d.z())*sqr(V.z());
}


Foam::symmTensor Foam::logStrainMechanicalConstitutiveLaw::expm
(
    const symmTensor& A
)
{
    tensor V(tensor::zero);
    vector d(vector::zero);
    eig3::eigen_decomposition(A, V, d);

    return
        Foam::exp(d.x())*sqr(V.x())
      + Foam::exp(d.y())*sqr(V.y())
      + Foam::exp(d.z())*sqr(V.z());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::logStrainMechanicalConstitutiveLaw::logStrainMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    subLawPtr_
    (
        mechanicalConstitutiveLaw::New(subLawDict(dict, "mechanicalLaw"))
    ),
    subStateName_("mechanicalLaw")
{
    if (!subLawPtr_->providesElasticPredictor())
    {
        FatalIOErrorInFunction(dict)
            << "The sub-law " << subLawPtr_->type() << " cannot be driven "
            << "by a trial elastic strain, so it cannot be used with "
            << type() << "." << nl
            << "Only an isotropic law with scalar-valued history can: it "
            << "must provide evaluateElasticPredictor."
            << exit(FatalIOError);
    }

    // Note: the planeStress entry is injected into this dictionary by the
    // mechanicalConstitutiveLawManager from the top-level entry in
    // mechanicalProperties
    if (Switch(dict.lookupOrDefault<Switch>("planeStress", false)))
    {
        FatalIOErrorInFunction(dict)
            << type() << " is not implemented for planeStress. If needed, "
            << "you can solve the case in 3-D and set the back to a "
            << "symmetryPlane and the front to traction-free"
            << exit(FatalIOError);
    }

    Info<< "    Logarithmic strain on top of " << subLawPtr_->type() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::logStrainMechanicalConstitutiveLaw::localConvergenceScale
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawState& state
) const
{
    // The elastic strain at the start of the step plus the strain of the
    // step, which bounds the trial elastic strain to first order. It only
    // scales a tolerance, so it is not worth the tensor logarithms the trial
    // strain itself would cost
    const UIndirectList<tensor>& F = kin.F();
    const UIndirectList<tensor>& Finv0 = kin.Finv0();

    const Field<symmTensor>& epsilonE0 = state.getSymmTensorField0("epsilonE");

    scalar maxMagEpsilonE = 0;

    forAll(F, i)
    {
        const tensor relF(F[i] & Finv0[i]);

        maxMagEpsilonE =
            max
            (
                maxMagEpsilonE,
                mag(epsilonE0[i]) + mag(symm(relF) - symmTensor(I))
            );
    }

    return maxMagEpsilonE;
}


const Foam::mechanicalConstitutiveLaw&
Foam::logStrainMechanicalConstitutiveLaw::childLaw
(
    const word& name
) const
{
    if (name != subStateName_)
    {
        FatalErrorInFunction
            << "Unknown child law '" << name << "'. This law owns only '"
            << subStateName_ << "'."
            << exit(FatalError);
    }

    return subLawPtr_();
}


void Foam::logStrainMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // The elastic left Cauchy-Green tensor, which carries the plastic
    // deformation from step to step: a restart without it would start again
    // from an unloaded state. Its undeformed value is the identity
    spec.addSymmTensor
    (
        "bE",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        symmTensor(I)
    );

    // The elastic logarithmic strain, 0.5*log(bE), kept for viewing and for
    // the convergence scale
    spec.addSymmTensor
    (
        "epsilonE",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        symmTensor::zero
    );
}


void Foam::logStrainMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    if (response.wantsVolumetricSplit())
    {
        FatalErrorInFunction
            << type() << " does not provide a volumetric split"
            << exit(FatalError);
    }

    UIndirectList<symmTensor>& sigma = response.stress();
    const UIndirectList<scalar>& J = kin.J();

    // Pre-processing: the trial elastic logarithmic strain, and its
    // eigenvectors
    Field<symmTensor> epsilonETrial(sigma.size());
    Field<tensor> eigenVectors(sigma.size());
    trialElasticStrain(kin, state, epsilonETrial, eigenVectors);

    // The sub-law's return mapping, from the trial strain, gives the elastic
    // strain the step ends at and the Kirchhoff stress, written straight
    // into the caller's storage
    Field<symmTensor>& epsilonE = state.symmTensorField("epsilonE");

    subLawPtr_->evaluateElasticPredictor
    (
        epsilonETrial,
        inputs,
        state.child(subStateName_),
        epsilonE,
        response
    );

    // Post-processing: the Cauchy stress, and the elastic left Cauchy-Green
    // tensor bE = exp(2*epsilonE). The return mapping of a law that
    // providesElasticPredictor() keeps the elastic strain coaxial with its
    // trial value, so the trial eigenvectors are its eigenvectors too, and
    // the exponential needs no decomposition of its own
    Field<symmTensor>& bE = state.symmTensorField("bE");

    forAll(sigma, i)
    {
        sigma[i] /= J[i];

        const tensor& V = eigenVectors[i];

        bE[i] =
            Foam::exp(2.0*(V.x() & epsilonE[i] & V.x()))*sqr(V.x())
          + Foam::exp(2.0*(V.y() & epsilonE[i] & V.y()))*sqr(V.y())
          + Foam::exp(2.0*(V.z() & epsilonE[i] & V.z()))*sqr(V.z());
    }

    // No analytical fourth-order tangent has been derived for this law
    fourthOrderByFiniteDifferenceOnly(kin, inputs, state, response);
}


// ************************************************************************* //
