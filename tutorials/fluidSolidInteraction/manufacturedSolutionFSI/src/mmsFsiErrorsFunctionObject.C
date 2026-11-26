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

#include "mmsFsiErrorsFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvcDdt.H"
#include "solidTractionFvPatchVectorField.H"
#include "mmsFsiFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mmsFsiErrorsFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mmsFsiErrorsFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Append the weighted L1 and L2 norms and the LInf norm of e to norms
static void appendNorms
(
    const scalarField& e,
    const scalarField& w,
    DynamicList<scalar>& norms
)
{
    const scalar sumW = gSum(w);

    norms.append(gSum(w*e)/sumW);
    norms.append(Foam::sqrt(gSum(w*sqr(e))/sumW));
    norms.append(gMax(e));
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::mmsFsiErrorsFunctionObject::writeErrors()
{
    const scalar t = time_.value();

    if (!filePtr_ && Pstream::master())
    {
        const fileName dir(time_.globalPath()/"postProcessing"/"mmsFsiErrors");
        mkDir(dir);
        filePtr_.reset(new OFstream(dir/"errors.dat"));

        filePtr_()
            << "# time"
            << " U_L1 U_L2 U_Linf p_L1 p_L2 p_Linf pMean_L2"
            << " D_L1 D_L2 D_Linf Us_L1 Us_L2 Us_Linf"
            << " intD_L1 intD_L2 intD_Linf"
            << " intTrac_L1 intTrac_L2 intTrac_Linf"
            << " intP_L1 intP_L2 intP_Linf" << endl;

        filePtr_().precision(10);
    }

    DynamicList<scalar> norms;

    // Fluid velocity and pressure at the current cell centres
    const fvMesh& fluidMesh = time_.lookupObject<fvMesh>(fluidRegion_);
    const volVectorField& U = fluidMesh.lookupObject<volVectorField>("U");
    const volScalarField& p = fluidMesh.lookupObject<volScalarField>("p");
    {
        const vectorField& C = fluidMesh.C().primitiveField();
        const scalarField& V = fluidMesh.V().field();

        scalarField eU(C.size());
        scalarField ep(C.size());
        forAll(C, cellI)
        {
            const mmsFsiFields::fluidState s =
                mmsFsiFields::fluid(C[cellI], t);

            eU[cellI] = mag(U[cellI] - s.U);
            ep[cellI] = p[cellI] - s.p;
        }

        appendNorms(eU, V, norms);
        appendNorms(mag(ep), V, norms);

        const scalar meanEp = gSum(V*ep)/gSum(V);
        norms.append(Foam::sqrt(gSum(V*sqr(ep - meanEp))/gSum(V)));

        if (writeFields_ && time_.writeTime())
        {
            volScalarField UError
            (
                IOobject
                (
                    "UError",
                    time_.timeName(),
                    fluidMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidMesh,
                dimensionedScalar("zero", U.dimensions(), 0)
            );
            UError.primitiveFieldRef() = eU;
            UError.write();

            volScalarField pError
            (
                IOobject
                (
                    "pError",
                    time_.timeName(),
                    fluidMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidMesh,
                dimensionedScalar("zero", p.dimensions(), 0)
            );
            pError.primitiveFieldRef() = ep;
            pError.write();
        }
    }

    // Solid displacement and velocity at the undeformed cell centres
    const fvMesh& solidMesh = time_.lookupObject<fvMesh>(solidRegion_);
    const volVectorField& D = solidMesh.lookupObject<volVectorField>("D");
    {
        const vectorField& C = solidMesh.C().primitiveField();
        const scalarField& V = solidMesh.V().field();
        const volVectorField Us(fvc::ddt(D));

        scalarField eD(C.size());
        scalarField eUs(C.size());
        forAll(C, cellI)
        {
            const mmsFsiFields::solidState s =
                mmsFsiFields::solid(C[cellI], t);

            eD[cellI] = mag(D[cellI] - s.D);
            eUs[cellI] = mag(Us[cellI] - s.U);
        }

        appendNorms(eD, V, norms);
        appendNorms(eUs, V, norms);
    }

    // Solid interface displacement and applied traction
    {
        const label patchID =
            solidMesh.boundaryMesh().findPatchID(interfacePatch_);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Patch " << interfacePatch_ << " not found in region "
                << solidRegion_ << exit(FatalError);
        }

        const fvPatchVectorField& patchD = D.boundaryField()[patchID];

        if (!isA<solidTractionFvPatchVectorField>(patchD))
        {
            FatalErrorInFunction
                << "The solid " << interfacePatch_ << " patch must be of "
                << "type solidTraction" << exit(FatalError);
        }

        const solidTractionFvPatchVectorField& tractionD =
            refCast<const solidTractionFvPatchVectorField>(patchD);

        const vectorField& Cf = patchD.patch().Cf();
        const scalarField& magSf = patchD.patch().magSf();

        scalarField eD(Cf.size());
        scalarField eTraction(Cf.size());
        forAll(Cf, faceI)
        {
            eD[faceI] =
                mag(patchD[faceI] - mmsFsiFields::solid(Cf[faceI], t).D);

            // The solid mesh is not moved, so the face centre gives the
            // reference abscissa; n is the fluid outward normal, i.e. minus
            // the solid outward normal
            vector x;
            vector tr;
            vector n;
            mmsFsiFields::interface(Cf[faceI].x(), t, x, tr, n);

            const vector appliedTraction =
                tractionD.traction()[faceI]
              + tractionD.pressure()[faceI]*n;

            eTraction[faceI] = mag(appliedTraction + tr);
        }

        appendNorms(eD, magSf, norms);
        appendNorms(eTraction, magSf, norms);
    }

    // Fluid interface pressure at the current face centres
    {
        const label patchID =
            fluidMesh.boundaryMesh().findPatchID(interfacePatch_);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Patch " << interfacePatch_ << " not found in region "
                << fluidRegion_ << exit(FatalError);
        }

        const fvPatchScalarField& patchP = p.boundaryField()[patchID];
        const vectorField& Cf = patchP.patch().Cf();

        scalarField ep(Cf.size());
        forAll(Cf, faceI)
        {
            ep[faceI] =
                mag(patchP[faceI] - mmsFsiFields::fluid(Cf[faceI], t).p);
        }

        appendNorms(ep, patchP.patch().magSf(), norms);
    }

    if (Pstream::master())
    {
        filePtr_() << t;
        forAll(norms, i)
        {
            filePtr_() << ' ' << norms[i];
        }
        filePtr_() << endl;
    }

    Info<< type() << ": U_L2 " << norms[1] << ", p_L2 " << norms[4]
        << ", D_L2 " << norms[8] << ", intTrac_L2 " << norms[17] << endl;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mmsFsiErrorsFunctionObject::mmsFsiErrorsFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    fluidRegion_("fluid"),
    solidRegion_("solid"),
    interfacePatch_("interface"),
    writeFields_(false),
    filePtr_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mmsFsiErrorsFunctionObject::read(const dictionary& dict)
{
    fluidRegion_ = dict.lookupOrDefault<word>("fluidRegion", "fluid");
    solidRegion_ = dict.lookupOrDefault<word>("solidRegion", "solid");
    interfacePatch_ =
        dict.lookupOrDefault<word>("interfacePatch", "interface");
    writeFields_ = dict.lookupOrDefault<Switch>("writeFields", false);

    return true;
}


bool Foam::mmsFsiErrorsFunctionObject::execute()
{
    return writeErrors();
}


bool Foam::mmsFsiErrorsFunctionObject::write()
{
    return true;
}


// ************************************************************************* //
