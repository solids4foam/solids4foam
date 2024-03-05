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

#include "abaqusUmatLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(abaqusUmatLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, abaqusUmatLinearElastic, linGeomMechLaw
    );

    // Declare fortran function prototypes
    extern "C"
    {
        // Note: all lowercase letters even if the fortran function has
        // uppercase letters
        void umat_
        (
            double STRESS[6],
            const double*, // double STATEV[6],
            double DDSDDE[6][6],
            const double*, // SSE,
            const double*, // SPD,
            const double*, // SCD,
            const double*, // RPL,
            const double*, // DDSDDT,
            const double*, // DRPLDE,
            const double*, // DRPLDT,
            const double STRAN[6],
            const double*, // const double DSTRAN[6],
            const double*, // TIME,
            const double*, // DTIME,
            const double*, // TEMP,
            const double*, // DTEMP,
            const double*, // PREDEF,
            const double*, // DPRED,
            const double*, // CMNAME,
            const int* NDI,
            const int* NSHR,
            const int* NTENS,
            const int* NSTATV,
            const double PROPS[],
            const int* NPROPS,
            const double*, // COORDS,
            const double*, // DROT,
            const double*, // PNEWDT,
            const double*, // CELENT,
            const double*, // DFGRD0,
            const double*, // DFGRD1,
            const double*, // NOEL,
            const double*, // NPT,
            const double*, // LAYER,
            const double*, // KSPT,
            const double*, // JSTEP,
            const double* // KINC
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::abaqusUmatLinearElastic::abaqusUmatLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    properties_(dict.lookup("properties")),
    // stateVariables_(0),
    impK_(dict.lookup("implicitStiffness")),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    )
{
    // Force storage of strain old time
    epsilon_.oldTime();

    // Initialise state varible fields

    // const scalarList stateVariablesInitialValues
    // (
    //     dict.lookup("stateVariablesInitialValues")
    // );

    // stateVariables_.setSize(stateVariablesInitialValues.size());

    // forAll(stateVariables_, fieldI)
    // {
    //     stateVariables_.set
    //     (
    //         fieldI,
    //         new volScalarField
    //         (
    //             IOobject
    //             (
    //                 "stateVariable" + Foam::name(fieldI + 1),
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::READ_IF_PRESENT,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh,
    //             dimensionedScalar
    //             (
    //                 "zero", dimless, stateVariablesInitialValues[fieldI]
    //             )
    //         )
    //     );

    //     // Force the old-time to be stored
    //     stateVariables_[fieldI].oldTime();
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::abaqusUmatLinearElastic::~abaqusUmatLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::abaqusUmatLinearElastic::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().correctBoundaryConditions();
#else
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


Foam::dimensionedScalar Foam::abaqusUmatLinearElastic::rhoScalar() const
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::abaqusUmatLinearElastic::impK() const
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
            impK_
        )
    );
}


void Foam::abaqusUmatLinearElastic::correct(volSymmTensorField& sigma)
{
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon_ = epsilon_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon_ = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::abaqusUmatLinearElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Not implemented for planeStress!" << abort(FatalError);
        }
    }

    // Initialise arrays to be passed to the fortran sub-routine, and then pass
    // them.
    // Note: many of the variables are not used by the abaqusUmatLinearElastic.C and
    // so we do not initialise them. In general, for other UMATS these may also
    // have to be initialised
    // BE CAREFUL: fortran expects column major indexing as opposed to row major
    // indexing, so we need to transpose all tensors before passing them and
    // receiving them
    // The following references are useful:
    // https://simplifiedfem.wordpress.com/about/tutorial-write-a-simple-umat-in-abaqus/
    // http://130.149.89.49:2080/v2016/books/sub/default.htm
    // http://130.149.89.49:2080/v2016/books/usb/default.htm?startat=pt01ch01s02aus02.html#usb-int-iconventions

    // Length of a stress tensor vector
    // Note: abaqus stores tensors as 1-D arrays (Voight notation_
    const int NTENS = 6;
    const int NDI = 3; // number of direct components
    const int NSHR = 3; // number of shear components
    const int NSTATV = 0; //stateVariablesI.size();

    // Material properties
    const int NPROPS = properties_.size();
    double PROPS[NPROPS];
    forAll(properties_, propI)
    {
        PROPS[propI] = properties_[propI];
    }

    // Internal field
    const symmTensorField& epsilonI = epsilon_.internalField();
#ifdef OPENFOAM_NOT_EXTEND
    symmTensorField& sigmaI = sigma.primitiveFieldRef();
#else
    symmTensorField& sigmaI = sigma.internalField();
#endif
    forAll(epsilonI, cellI)
    {
        double STRESS[NTENS];
        // Initialise stress to zero
        for (int i = 0; i < NTENS; i++)
        {
            STRESS[i] = 0.0;
        }
        // double STATEV[NSTATV];
        // forAll(stateVariables_, varI)
        // {
        //     STATEV[varI] = stateVariables_.internalField()[varI];
        // }
        double DDSDDE[NTENS][NTENS];
        // double SSE[];
        // double SPD[];
        // double SCD[];
        // double RPL[];
        // double DDSDDT[NTENS];
        // double DRPLDE[NTENS];
        // double DRPLDT[];
        const symmTensor& eps = epsilonI[cellI];
        double STRAN[NTENS] =
        {
            eps.xx(), eps.yy(), eps.zz(),
            eps.xy(), eps.yz(), eps.xz()
        };
        // double DSTRAN[NTENS];
        // double TIME[2];
        // double DTIME;
        // double TEMP[];
        // double DTEMP[];
        // double PREDEF;
        // double DPRED;
        // double CMNAME[];
        // double COORDS[];
        // double DROT[];
        // double PNEWDT[];
        // double CELENT[];
        // double DFGRD0[];
        // double DFGRD1[];
        // double NOEL[];
        // double NPT[];
        // double LAYER[];
        // double KSPT[];
        // double JSTEP[4];
        // double KINC[];

        // Call Abaqus UMAT to calculate the stress
        double notImplemented = 0;
        umat_
        (
            STRESS,
            &notImplemented, // STATEV,
            DDSDDE,
            &notImplemented, // SSE,
            &notImplemented, // SPD,
            &notImplemented, // SCD,
            &notImplemented, // RPL,
            &notImplemented, // DDSDDT,
            &notImplemented, // DRPLDE,
            &notImplemented, // DRPLDT,
            STRAN,
            &notImplemented, //DSTRAN,
            &notImplemented, // TIME,
            &notImplemented, // DTIME,
            &notImplemented, // TEMP,
            &notImplemented, // DTEMP,
            &notImplemented, // PREDEF,
            &notImplemented, // DPRED,
            &notImplemented, // CMNAME,
            &NDI,
            &NSHR,
            &NTENS,
            &NSTATV,
            PROPS,
            &NPROPS,
            &notImplemented, // COORDS,
            &notImplemented, // DROT,
            &notImplemented, // PNEWDT,
            &notImplemented, // CELENT,
            &notImplemented, // DFGRD0,
            &notImplemented, // DFGRD1,
            &notImplemented, // NOEL,
            &notImplemented, // NPT,
            &notImplemented, // LAYER,
            &notImplemented, // KSPT,
            &notImplemented, // JSTEP,
            &notImplemented // KINC
        );

        // Retrieve the stress
        // If used, you would also retrieve the state variables
        sigmaI[cellI].xx() = STRESS[0];
        sigmaI[cellI].yy() = STRESS[1];
        sigmaI[cellI].zz() = STRESS[2];
        sigmaI[cellI].xy() = STRESS[3];
        sigmaI[cellI].yz() = STRESS[4];
        sigmaI[cellI].xz() = STRESS[5];
    }

    // Boundary field
    forAll(epsilon_.boundaryField(), patchI)
    {
        const symmTensorField& epsilonP = epsilon_.boundaryField()[patchI];
#ifdef OPENFOAM_NOT_EXTEND
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
#else
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
#endif
        forAll(epsilonP, faceI)
        {
            double STRESS[NTENS];
            // Initialise stress to zero
            for (int i = 0; i < NTENS; i++)
            {
                STRESS[i] = 0.0;
            }
            // double STATEV[NSTATV];
            // forAll(stateVariables_, varI)
            // {
            //     STATEV[varI] = stateVariables_.internalField()[varI];
            // }
            double DDSDDE[NTENS][NTENS];
            // double SSE[];
            // double SPD[];
            // double SCD[];
            // double RPL[];
            // double DDSDDT[NTENS];
            // double DRPLDE[NTENS];
            // double DRPLDT[];
            const symmTensor& eps = epsilonP[faceI];
            double STRAN[NTENS] =
                {
                    eps.xx(), eps.yy(), eps.zz(),
                    eps.xy(), eps.yz(), eps.xz()
                };
            // double DSTRAN[NTENS];
            // double TIME[2];
            // double DTIME;
            // double TEMP[];
            // double DTEMP[];
            // double PREDEF;
            // double DPRED;
            // double CMNAME[];
            // double COORDS[];
            // double DROT[];
            // double PNEWDT[];
            // double CELENT[];
            // double DFGRD0[];
            // double DFGRD1[];
            // double NOEL[];
            // double NPT[];
            // double LAYER[];
            // double KSPT[];
            // double JSTEP[4];
            // double KINC[];

            // Call Abaqus UMAT to calculate the stress
            double notImplemented = 0;
            umat_
            (
                STRESS,
                &notImplemented, // STATEV,
                DDSDDE,
                &notImplemented, // SSE,
                &notImplemented, // SPD,
                &notImplemented, // SCD,
                &notImplemented, // RPL,
                &notImplemented, // DDSDDT,
                &notImplemented, // DRPLDE,
                &notImplemented, // DRPLDT,
                STRAN,
                &notImplemented, //DSTRAN,
                &notImplemented, // TIME,
                &notImplemented, // DTIME,
                &notImplemented, // TEMP,
                &notImplemented, // DTEMP,
                &notImplemented, // PREDEF,
                &notImplemented, // DPRED,
                &notImplemented, // CMNAME,
                &NDI,
                &NSHR,
                &NTENS,
                &NSTATV,
                PROPS,
                &NPROPS,
                &notImplemented, // COORDS,
                &notImplemented, // DROT,
                &notImplemented, // PNEWDT,
                &notImplemented, // CELENT,
                &notImplemented, // DFGRD0,
                &notImplemented, // DFGRD1,
                &notImplemented, // NOEL,
                &notImplemented, // NPT,
                &notImplemented, // LAYER,
                &notImplemented, // KSPT,
                &notImplemented, // JSTEP,
                &notImplemented // KINC
            );

            // Retrieve the stress
            // If used, you would also retrieve the state variables
            sigmaP[faceI].xx() = STRESS[0];
            sigmaP[faceI].yy() = STRESS[1];
            sigmaP[faceI].zz() = STRESS[2];
            sigmaP[faceI].xy() = STRESS[3];
            sigmaP[faceI].yz() = STRESS[4];
            sigmaP[faceI].xz() = STRESS[5];
        }
    }
}

void Foam::abaqusUmatLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    notImplemented("Foam::abaqusUmatLinearElastic::correct(surfaceSymmTensorField)");
}


// ************************************************************************* //
