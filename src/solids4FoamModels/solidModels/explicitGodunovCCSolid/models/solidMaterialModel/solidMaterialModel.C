/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidMaterialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidMaterialModel, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

solidMaterialModel::solidMaterialModel
(
    const volTensorField& F,
    const dictionary& dict
)
:
    P_
    (
        IOobject
        (
            "P",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedTensor("P", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
    ),

    p_
    (
        IOobject
        (
            "p",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    energyAlgorithm_
    (
        IOobject
        (
            "energyAlgorithm",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedScalar
        (
            "energyAlgorithm",
            dimensionSet(1,-1,-2,0,0,0,0),
            0.0
        )
    ),

    rho_("rho",dimDensity , 0.0),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    mu_("mu", dimPressure, 0.0),
    lambda_("lambda", dimPressure, 0.0),
    kappa_("kappa_", dimPressure, 0.0),
    Up_("Up_", dimVelocity, 0.0),
    Us_("Us_", dimVelocity, 0.0)
{
    // Read the mechanical laws
    const PtrList<entry> lawEntries(dict.lookup("mechanical"));

    const dictionary& materialDict = lawEntries[0].dict();

    // Read model rho, E, and nu from the material dictionary

    const word model
    (
        materialDict.lookup("type")
    );
    model_ = model;

    rho_ = dimensionedScalar(materialDict.lookup("rho"));
    E_ = dimensionedScalar(materialDict.lookup("E"));
    nu_ = dimensionedScalar(materialDict.lookup("nu"));
    mu_ =(E_/(2.0*(1.0 + nu_)));

    lambda_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_)));

    kappa_ =lambda_ + (2.0/3.0)*mu_;

    Up_ = sqrt((lambda_+2.0*mu_)/rho_);
    Us_ = sqrt(mu_/rho_);


    correct();// add here form creatField.H in original solver to make sure variables are apdated like rho_ ...
    p_.write();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
solidMaterialModel::~solidMaterialModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidMaterialModel::correct()
{
    const fvMesh& mesh_ = P_.mesh();
    const objectRegistry& db = mesh_.thisDb();
    const volTensorField& H_ = db.lookupObject<volTensorField>("H");
    const volTensorField& F_ = db.lookupObject<volTensorField>("F");
    const volScalarField& J_ = db.lookupObject<volScalarField>("J");
    const volVectorField& lm_ = db.lookupObject<volVectorField>("lm");

    energyAlgorithm_ =
        0.5*mu_*(pow(J_,(-2.0/3.0))*(F_ && F_)-3.0)
      + (0.5*kappa_*(J_ - 1.0)*(J_ - 1.0)) + (0.5*(lm_ & lm_)/rho_);

    if (model_ == "linearElastic")
    {
        p_ = kappa_*(tr(F_) - 3.0);
        P_ = mu_*(F_ + F_.T() - ((2.0/3.0)*tr(F_)*tensor::I)) + p_*tensor::I;
    }

    else if (model_ == "neoHookeanElastic")
    {
        p_ = kappa_*(J_-1.0);
        P_ =
            mu_*pow(J_,(-2.0/3.0))*F_
          - ( (mu_/3.0)*pow(J_,(-5.0/3.0))*(F_ && F_)*H_ ) + p_*H_;
    }

    else
    {
        FatalErrorIn
        (
            "solidMaterialModel.C"
        )   << "Valid type entries are 'linearElastic' or 'neoHookeanElastic' for"
            << "solidMaterialModel"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void solidMaterialModel::printMaterialProperties()
{
    Info<< "\nPrinting material properties ..." << nl
        << "Constitutive model = " << model_ << nl
        << "Density = " << rho_.value() << " " << rho_.dimensions() << nl
        << "Young's modulus = " << E_.value() << " " << E_.dimensions() << nl
        << "Poisson's ratio = " << nu_.value() << " " << nu_.dimensions() << nl
        << "Lame's first parameter lambda = " << lambda_.value() << " "
        << lambda_.dimensions() << nl
        << "Lame's second parameter mu = " << mu_.value() << " "
        << mu_.dimensions() << nl
        << "Bulk modulus kappa = " << kappa_.value() << " "
        << kappa_.dimensions() << nl
        << "Linear pressure wave speed = " << Up_.value() << " "
        << Up_.dimensions() << nl
        << "Linear shear wave speed = " << Us_.value() << " "
        << Us_.dimensions() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //