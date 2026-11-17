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

Application
    stabilisationFourierCheck

Description
    Verifies the runtime-selected pressure-stabilisation models against
    independent analytical Fourier symbols (fourierSymbols.H) on a uniform,
    axis-aligned periodic Cartesian mesh.

    Each model is applied, through the calls a solid model makes, to discrete
    Fourier modes of the pressure. The measured eigenvalue (Rayleigh quotient)
    and the eigenmode residual are compared with the analytical symbol, with
    and without spectral normalisation. Results are written as CSV files to
    postProcessing/stabilisationFourierCheck, and the utility exits nonzero
    if a check fails.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "stabilisationModel.H"
#include "compatibilityFunctions.H"
#include "cartesianMeshInfo.H"
#include "fourierSymbols.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "mathematicalConstants.H"

// std::isfinite and std::lround: no OpenFOAM scalar equivalent exists in all
// supported flavours. std::memcmp: the bitwise-equality checks must
// distinguish -0 from +0 and compare NaN payloads, which == does not
#include <cmath>
#include <cstring>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static bool optionFound(const argList& args, const word& option)
{
#ifdef OPENFOAM_COM
    return args.found(option);
#else
    return args.optionFound(option);
#endif
}


template<class Type>
static void readOption(const argList& args, const word& option, Type& value)
{
#ifdef OPENFOAM_COM
    args.readIfPresent(option, value);
#else
    args.optionReadIfPresent(option, value);
#endif
}


static scalar piValue()
{
#ifdef FOAMEXTEND
    return mathematicalConstant::pi;
#else
    return constant::mathematical::pi;
#endif
}


static label modelIndex(const wordList& names, const word& name)
{
    forAll(names, i)
    {
        if (names[i] == name)
        {
            return i;
        }
    }

    return -1;
}


//- True for the models whose spectral shape is directionally split
static bool splitSpectralShape(const word& type)
{
    return type == "diffStencilLaplacian" || type == "RhieChow";
}


//- True for a model that applies the plain Laplacian: laplacian, or
//  generalisedEvenOrderLaplacian with laplacianPower 0
static bool plainLaplacian(const dictionary& modelDict)
{
    const word type(modelDict.lookup("type"));

    return
        type == "laplacian"
     || (
            type == "generalisedEvenOrderLaplacian"
         && readLabel(modelDict.lookup("laplacianPower")) == 0
        );
}


//- True if two models share their implementation, so that their eigenvalues
//  must agree bitwise: two plain Laplacians, or two directionally split
//  models. The models are identified by their entries, not their names, so
//  the order of a pair does not matter
static bool sharedImplementation(const dictionary& a, const dictionary& b)
{
    return
        (plainLaplacian(a) && plainLaplacian(b))
     || (
            splitSpectralShape(word(a.lookup("type")))
         && splitSpectralShape(word(b.lookup("type")))
        );
}


//- Number of repeated Laplacian factors, m + 1 with m = laplacianPower, or
//  -1 for the directionally split models
static label nLaplacianFactors(const word& type, const dictionary& modelDict)
{
    if (splitSpectralShape(type))
    {
        return -1;
    }

    if (type == "JamesonSchmidtTurkel")
    {
        return 2;
    }

    if (type == "generalisedEvenOrderLaplacian")
    {
        return readLabel(modelDict.lookup("laplacianPower")) + 1;
    }

    return 1;
}


//- Name of the symbol family that the analytical eigenvalue belongs to
static word spectralRelation
(
    const word& type,
    const dictionary& modelDict,
    const bool normalise
)
{
    if (splitSpectralShape(type))
    {
        return "splitM2Shape";
    }

    if
    (
        !normalise
     && (
            type == "laplacian"
         || (
                type == "generalisedEvenOrderLaplacian"
             && readLabel(modelDict.lookup("laplacianPower")) == 0
            )
        )
    )
    {
        return "legacyM1";
    }

    return "coupledRepeatedLaplacian";
}


//- Record the model settings, including nested sub-dictionaries, as
//  key=value; pairs
static string parameters(const dictionary& dict, const string& prefix = "")
{
    string result;
    const wordList keys(dict.toc());

    forAll(keys, i)
    {
        const string key(prefix + keys[i]);

        if (dict.isDict(keys[i]))
        {
            result += parameters(dict.subDict(keys[i]), key + ".");
        }
        else
        {
            const ITstream& tokens = dict.lookup(keys[i]);
            OStringStream value;

            forAll(tokens, j)
            {
                if (j)
                {
                    value << ' ';
                }
                value << tokens[j];
            }

            result += key + "=" + value.str().c_str() + ";";
        }
    }

    return result;
}


//- Quote a CSV field, doubling any embedded quote
static string csvQuote(const string& value)
{
    string result("\"");

    for (const char c : value)
    {
        if (c == '"')
        {
            result += '"';
        }
        result += c;
    }

    return result + "\"";
}


//- Open a CSV file and write the common metadata lines. Strings are written
//  as const char*: Ostream would quote a string or fileName
static autoPtr<OFstream> openCsv
(
    const fileName& path,
    const Time& runTime,
    const scalar tolerance,
    const scalar pairTolerance,
    const scalar gammaTolerance
)
{
    autoPtr<OFstream> osPtr(new OFstream(path));
    OFstream& os = osPtr();

    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open " << path << " for writing" << exit(FatalError);
    }

    string version(getEnv("WM_PROJECT_VERSION"));
    if (version.empty())
    {
        version = "unknown";
    }

    os.precision(17);
    os  << "# OpenFOAM=" << version.c_str()
        << "\n# case=" << runTime.path().c_str()
        << "\n# gamma=1 (rAUf dimensions area/pressure)"
        << "\n# tolerance=" << tolerance
        << "\n# pairTolerance=" << pairTolerance
        << "\n# gammaTolerance=" << gammaTolerance << '\n';

    return osPtr;
}


//- Flush a CSV file and return whether everything was written
static bool closeCsv(OFstream& os)
{
    os.flush();

    return os.good();
}


//- Validate a mode, given as theta/pi, against the lattice of the mesh and
//  return its angles theta, with the integer wave numbers in k
static vector fourierModeAngles
(
    const fvMesh& mesh,
    const cartesianMeshInfo& grid,
    const vector& thetaByPi,
    const string& modeName,
    label k[3]
)
{
    const scalar pi = piValue();
    vector theta(vector::zero);

    for (direction d = 0; d < 3; ++d)
    {
        // The wave number k_d = N_d*(theta_d/pi)/2 must be an integer in
        // [0, N_d/2], and zero in the empty direction
        const scalar discrete = grid.N[d]*thetaByPi[d]/2;

        if
        (
            !std::isfinite(discrete)
         || discrete < 0
         || discrete > grid.N[d]/2.0
        )
        {
            FatalErrorInFunction
                << modeName.c_str() << " must lie in [0,pi]"
                << exit(FatalError);
        }

        k[d] = std::lround(discrete);

        if
        (
            mag(discrete - k[d]) > 1e-12
         || (mesh.geometricD()[d] == -1 && thetaByPi[d] != 0)
        )
        {
            FatalErrorInFunction
                << modeName.c_str() << " is incompatible with mesh periodicity"
                << exit(FatalError);
        }

        theta[d] = 2*pi*k[d]/grid.N[d];
    }

    return theta;
}


//- Set p to the Fourier mode cos(theta & i) of the cell indices i. Using the
//  index, rather than the cell centre, keeps the Nyquist mode nonzero:
//  cos(pi*(i + 1/2)) would vanish in every cell
static void setFourierPressure
(
    volScalarField& p,
    const cartesianMeshInfo& grid,
    const vector& theta
)
{
    scalarField& pI = primitiveFieldRef(p);

    forAll(pI, celli)
    {
        scalar phase = 0;
        for (direction d = 0; d < 3; ++d)
        {
            phase += theta[d]*grid.index[d][celli];
        }

        pI[celli] = Foam::cos(phase);
    }

    p.correctBoundaryConditions();
}


//- Construct each configured model and return; used by the rejection tests
static void constructModels
(
    const fvMesh& mesh,
    const dictionary& models,
    const wordList& names
)
{
    forAll(names, modeli)
    {
        autoPtr<stabilisationModel> model
        (
            stabilisationModel::New
            (
                mesh,
                models.subDict(names[modeli]),
                dimPressure/dimLength
            )
        );

        Info<< "Constructed stabilisation model " << names[modeli] << nl;
    }
}


//- Generic field export (-writeFields with a fieldOutput entry): write the
//  pressure and stabilisation fields of each named mode and model. Returns
//  the number of modes that are not eigenmodes to within the tolerance
static label writeFieldOutput
(
    const fvMesh& mesh,
    const cartesianMeshInfo& grid,
    const dictionary& dict,
    const dictionary& models,
    const wordList& names,
    const scalar tol
)
{
    const Time& runTime = mesh.time();
    const dictionary& outputDict = dict.subDict("fieldOutput");

    if (!outputDict.found("modes"))
    {
        FatalIOErrorInFunction(dict)
            << "fieldOutput requires a modes sub-dictionary"
            << exit(FatalIOError);
    }

    const dictionary& outputModes = outputDict.subDict("modes");
    const wordList modeNames(outputModes.toc());

    if (names.empty() || modeNames.empty())
    {
        FatalIOErrorInFunction(dict)
            << "fieldOutput models and modes must not be empty"
            << exit(FatalIOError);
    }

    const word pressurePrefix
    (
        outputDict.lookupOrDefault<word>("pressurePrefix", "fourierP")
    );
    const word stabilisationPrefix
    (
        outputDict.lookupOrDefault<word>("stabilisationPrefix", "fourierS")
    );

    wordList patchTypes(mesh.boundary().size());
    forAll(patchTypes, patchi)
    {
        patchTypes[patchi] = mesh.boundary()[patchi].type();
    }

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0),
        patchTypes
    );

    surfaceScalarField rAUf
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimArea/dimPressure, 1)
    );

    Info<< "Fourier field output: "
        << grid.N[0] << 'x' << grid.N[1] << 'x' << grid.N[2]
        << ", " << modeNames.size() << " modes, "
        << names.size() << " models" << nl;

    label failures = 0;

    forAll(modeNames, modei)
    {
        const word& modeName = modeNames[modei];
        label k[3];
        const vector theta
        (
            fourierModeAngles
            (
                mesh,
                grid,
                vector(outputModes.lookup(modeName)),
                "Field-output mode " + modeName,
                k
            )
        );

        setFourierPressure(p, grid, theta);

        volScalarField pMode
        (
            IOobject
            (
                pressurePrefix + "_" + modeName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p
        );
        pMode.write();

        volVectorField gradP(fvc::grad(p));
        gradP.correctBoundaryConditions();

        const scalarField pressure(primitiveField(p));
        const scalar p2 = gSum(mesh.V()*sqr(pressure));

        if (p2 <= SMALL)
        {
            FatalErrorInFunction
                << "Degenerate field-output mode " << modeName
                << exit(FatalError);
        }

        forAll(names, modeli)
        {
            autoPtr<stabilisationModel> model
            (
                stabilisationModel::New
                (
                    mesh,
                    models.subDict(names[modeli]),
                    dimPressure/dimLength
                )
            );

            model->updateScalar(p, &gradP);

            const volScalarField& stab = model->cellScalar(&rAUf, true);
            const scalarField& S = primitiveField(stab);

            // Rayleigh quotient and the residual of the eigenmode relation
            const scalar num = gSum(mesh.V()*pressure*S)/p2;
            const scalar residual =
                Foam::sqrt(gSum(mesh.V()*sqr(S - num*pressure)));
            const scalar modeErr =
                residual/max(mag(num)*Foam::sqrt(p2), SMALL);

            failures += !std::isfinite(num) || modeErr > tol;

            volScalarField outputS
            (
                IOobject
                (
                    stabilisationPrefix + "_" + names[modeli] + "_" + modeName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                stab
            );
            outputS.write();

            Info<< modeName << ' ' << names[modeli]
                << " lambda_num=" << num
                << " eigenmodeError=" << modeErr << nl;
        }
    }

    Info<< "Fourier fields written to "
        << runTime.path()/runTime.timeName() << nl;

    return failures;
}


//- Write the analytical symbols along three cuts of the (thetaX, thetaY)
//  plane and on a 2-D map, for plotting. Returns false on a write error
static bool writeSymbolExport
(
    const fileName& outDir,
    const word& resolution,
    const Time& runTime,
    const cartesianMeshInfo& grid,
    const dictionary& models,
    const wordList& names,
    const scalar tol,
    const scalar pairTol,
    const scalar gammaTol
)
{
    const scalar pi = piValue();

    autoPtr<OFstream> sweepPtr
    (
        openCsv
        (
            outDir/("symbols_" + resolution + ".csv"),
            runTime,
            tol,
            pairTol,
            gammaTol
        )
    );
    autoPtr<OFstream> mapsPtr
    (
        openCsv
        (
            outDir/("symbolMaps_" + resolution + ".csv"),
            runTime,
            tol,
            pairTol,
            gammaTol
        )
    );
    OFstream& sweep = sweepPtr();
    OFstream& maps = mapsPtr();

    sweep
        << "model,cut,theta_over_pi,lambda_exact,hx,hy,scaleFactor\n";
    maps
        << "model,thetaX_over_pi,thetaY_over_pi,lambda_exact,hx,hy,"
        << "scaleFactor\n";

    forAll(names, modeli)
    {
        const dictionary& md = models.subDict(names[modeli]);
        const scalar scale = md.lookupOrDefault<scalar>("scaleFactor", 1);

        // Cuts: 0 along thetaX (thetaY = 0), 1 along the diagonal, 2 along
        // thetaY at thetaX = pi
        for (label cut = 0; cut < 3; ++cut)
        {
            for (label i = 0; i <= 200; ++i)
            {
                const scalar t = pi*i/200;
                const vector theta(cut == 2 ? pi : t, cut == 0 ? 0 : t, 0);

                sweep
                    << names[modeli] << ',' << cut << ',' << t/pi << ','
                    << expectedEigenvalue(md, theta, grid.h, grid.N) << ','
                    << grid.h.x() << ',' << grid.h.y() << ',' << scale
                    << '\n';
            }
        }

        for (label j = 0; j <= 64; ++j)
        {
            for (label i = 0; i <= 64; ++i)
            {
                const vector theta(pi*i/64, pi*j/64, 0);

                maps
                    << names[modeli] << ',' << i/64.0 << ',' << j/64.0 << ','
                    << expectedEigenvalue(md, theta, grid.h, grid.N) << ','
                    << grid.h.x() << ',' << grid.h.y() << ',' << scale
                    << '\n';
            }
        }
    }

    const bool sweepWritten = closeCsv(sweep);
    const bool mapsWritten = closeCsv(maps);

    return sweepWritten && mapsWritten;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#ifdef FOAMEXTEND
    argList::validOptions.insert("dict", "name");
    argList::validOptions.insert("writeFields", "");
    argList::validOptions.insert("symbolSweep", "");
    argList::validOptions.insert("gammaLinearityCheck", "scalar");
    argList::validOptions.insert("constructOnly", "");
#else
    argList::addOption
    (
        "dict",
        "name",
        "Dictionary in system (default stabilisationFourierCheckDict)"
    );
    argList::addBoolOption
    (
        "writeFields",
        "Write Fourier pressure and production stabilisation fields"
    );
    argList::addBoolOption
    (
        "symbolSweep",
        "Export analytical cuts and 2-D maps"
    );
    argList::addOption
    (
        "gammaLinearityCheck",
        "scalar",
        "Override secondary constant gamma (0 disables)"
    );
    argList::addBoolOption
    (
        "constructOnly",
        "Construct configured models, then exit"
    );
#endif

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Settings

    word dictName("stabilisationFourierCheckDict");
    readOption(args, "dict", dictName);

    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar tol = dict.lookupOrDefault<scalar>("tolerance", 1e-9);
    const scalar pairTol =
        dict.lookupOrDefault<scalar>("pairTolerance", 1e-12);
    const scalar gammaTol =
        dict.lookupOrDefault<scalar>("gammaTolerance", 1e-14);
    const Switch normalisedVerification
    (
        dict.lookupOrDefault<Switch>("normalisedVerification", false)
    );
    const label rStar = normalisedVerification
      ? readLabel(dict.lookup("referenceNyquistDirections"))
      : -1;

    scalar gammaCheck =
        dict.lookupOrDefault<scalar>("gammaLinearityCheck", 0.37);
    readOption(args, "gammaLinearityCheck", gammaCheck);

    if
    (
        !std::isfinite(tol) || tol <= 0
     || !std::isfinite(pairTol) || pairTol <= 0
     || !std::isfinite(gammaTol) || gammaTol <= 0
     || !std::isfinite(gammaCheck) || gammaCheck < 0
    )
    {
        FatalErrorInFunction
            << "Invalid tolerance or constant gamma"
            << exit(FatalError);
    }

    // Validates the mesh: uniform, axis-aligned, periodic
    const cartesianMeshInfo grid(mesh);
    const scalar pi = piValue();
    const dictionary& models = dict.subDict("models");
    const wordList names(models.toc());

    if (optionFound(args, "constructOnly"))
    {
        constructModels(mesh, models, names);

        Info<< "End" << endl;
        return 0;
    }

    if (optionFound(args, "writeFields") && dict.found("fieldOutput"))
    {
        const label failures =
            writeFieldOutput(mesh, grid, dict, models, names, tol);

        Info<< "Failures: " << failures << nl << "End" << endl;
        return failures ? 1 : 0;
    }

    const List<vector> modes(dict.lookup("modes"));
    const List<Pair<word>> pairs(dict.lookup("crossChecks"));

    if (names.empty() || modes.empty())
    {
        FatalErrorInFunction
            << "Models and modes must not be empty"
            << exit(FatalError);
    }

    // Output files

    const word resolution
    (
        Foam::name(grid.N[0]) + "x" + Foam::name(grid.N[1])
      + "x" + Foam::name(grid.N[2])
    );
    const fileName outDir
    (
        runTime.path()/"postProcessing"/"stabilisationFourierCheck"
    );
    mkDir(outDir);

    autoPtr<OFstream> resultsPtr
    (
        openCsv
        (
            outDir/("results_" + resolution + ".csv"),
            runTime,
            tol,
            pairTol,
            gammaTol
        )
    );
    autoPtr<OFstream> crossPtr
    (
        openCsv
        (
            outDir/("crossChecks_" + resolution + ".csv"),
            runTime,
            tol,
            pairTol,
            gammaTol
        )
    );
    autoPtr<OFstream> gammaPtr
    (
        openCsv
        (
            outDir/("gammaChecks_" + resolution + ".csv"),
            runTime,
            tol,
            pairTol,
            gammaTol
        )
    );
    OFstream& results = resultsPtr();
    OFstream& cross = crossPtr();
    OFstream& gamma = gammaPtr();

    // The paperOrder column is the number of repeated Laplacian factors,
    // m + 1 with m = laplacianPower (zero-based); the name is kept for the
    // existing post-processing of these files
    results
        << "spectralNormalisation,referenceNyquistDirections,model,type,"
        << "spectralRelation,paperOrder,scaleFactor,laplacianPower,params,"
        << "Nx,Ny,Nz,hx,hy,hz,modeIndex,kx,ky,kz,"
        << "thetaX_over_pi,thetaY_over_pi,thetaZ_over_pi,"
        << "lambda_exact,lambda_num,lambda_ref,"
        << "relativeEigenvalueError,rayleighErr,eigenmodeErr,"
        << "eigenmodeErrRef,referenceMode,pass\n";
    cross
        << "spectralNormalisation,referenceNyquistDirections,modelA,modelB,"
        << "modeIndex,thetaX_over_pi,thetaY_over_pi,maxAbsDiffOverRef,"
        << "bitwiseEqual,pass\n";
    gamma
        << "spectralNormalisation,referenceNyquistDirections,model,"
        << "modeIndex,gamma,maxAbsDiffOverRef,nullGammaBitwiseEqual,pass\n";

    // Fields: the pressure mode, and face weights gamma = rAUf = 1 and a
    // second constant for the linearity check

    wordList patchTypes(mesh.boundary().size());
    forAll(patchTypes, patchi)
    {
        patchTypes[patchi] = mesh.boundary()[patchi].type();
    }

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0),
        patchTypes
    );

    surfaceScalarField rAUf
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimArea/dimPressure, 1)
    );

    surfaceScalarField gammaC
    (
        IOobject
        (
            "gammaC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("c", dimArea/dimPressure, gammaCheck)
    );

    const scalarField& V = mesh.V();
    label failures = 0;

    // Campaign 0 checks the models as configured; campaign 1, with
    // normalisedVerification, checks them again with spectral normalisation
    // switched on
    const label nCampaigns = normalisedVerification ? 2 : 1;

    for (label campaigni = 0; campaigni < nCampaigns; ++campaigni)
    {
        const bool normalise = campaigni == 1;

        // Reference mode: Nyquist in every active direction for the raw
        // campaign, in the first rStar active directions when normalised
        vector referenceTheta(vector::zero);
        label referenceDirections = 0;

        for (direction d = 0; d < 3; ++d)
        {
            if
            (
                mesh.geometricD()[d] != -1
             && (!normalise || referenceDirections < rStar)
            )
            {
                referenceTheta[d] = pi;
                ++referenceDirections;
            }
        }

        if (normalise && referenceDirections != rStar)
        {
            FatalErrorInFunction
                << "referenceNyquistDirections exceeds active mesh directions"
                << exit(FatalError);
        }

        Info<< (normalise ? "NORMALISED" : "RAW") << " Fourier check: "
            << resolution << ", h=" << grid.h << ", tolerance=" << tol;
        if (normalise)
        {
            Info<< ", referenceNyquistDirections=" << rStar;
        }
        Info<< nl
            << "model mode(kx,ky,kz) theta/pi lambda_exact lambda_num "
            << "relativeEigenvalueError eigenmodeErr status" << nl;

        forAll(modes, modei)
        {
            OStringStream modeName;
            modeName << "Mode " << modes[modei];

            label k[3];
            const vector theta
            (
                fourierModeAngles(mesh, grid, modes[modei], modeName.str(), k)
            );

            setFourierPressure(p, grid, theta);

            volVectorField gradP(fvc::grad(p));
            gradP.correctBoundaryConditions();

            const scalarField pressure(primitiveField(p));
            const scalar p2 = gSum(V*sqr(pressure));

            if (p2 <= SMALL)
            {
                FatalErrorInFunction
                    << "Degenerate Fourier mode" << exit(FatalError);
            }

            // Measured responses and reference eigenvalue magnitudes of all
            // models, kept for the pair checks below
            List<scalarField> numerical(names.size());
            scalarField refs(names.size());

            forAll(names, modeli)
            {
                dictionary modelDict(models.subDict(names[modeli]));

                if (normalise)
                {
                    modelDict.set("spectralNormalisation", Switch(true));
                    modelDict.set("referenceNyquistDirections", rStar);
                }

                const word type(modelDict.lookup("type"));
                const scalar scale = type == "laplacian"
                  ? modelDict.lookupOrDefault<scalar>("scaleFactor", 1)
                  : readScalar(modelDict.lookup("scaleFactor"));
                const label power = type == "generalisedEvenOrderLaplacian"
                  ? readLabel(modelDict.lookup("laplacianPower"))
                  : -1;

                // Analytical eigenvalue of this mode, and the magnitude of
                // that of the reference mode, which scales the errors
                const scalar exact =
                    expectedEigenvalue(modelDict, theta, grid.h, grid.N);
                const scalar ref = mag
                (
                    expectedEigenvalue
                    (
                        modelDict,
                        referenceTheta,
                        grid.h,
                        grid.N
                    )
                );

                if (!std::isfinite(ref) || ref <= 0)
                {
                    FatalErrorInFunction
                        << "Verification needs a finite nonzero reference "
                        << "response" << exit(FatalError);
                }
                refs[modeli] = ref;

                // Measured response: the runtime-selected model, through the
                // updateScalar and cellScalar calls a solid model makes
                autoPtr<stabilisationModel> model
                (
                    stabilisationModel::New
                    (
                        mesh,
                        modelDict,
                        dimPressure/dimLength
                    )
                );
                model->updateScalar(p, &gradP);

                const volScalarField& stab = model->cellScalar(&rAUf, true);
                numerical[modeli] = primitiveField(stab);
                const scalarField& S = numerical[modeli];

                // Rayleigh quotient, and the residual of S = exact*p, which
                // is zero only if p is an eigenmode with that eigenvalue
                const scalar num = gSum(V*pressure*S)/p2;
                const scalar residual =
                    Foam::sqrt(gSum(V*sqr(S - exact*pressure)));
                const scalar rayleighErr = mag(num - exact)/ref;
                const scalar modeErrRef = residual/(ref*Foam::sqrt(p2));
                const scalar relative =
                    exact != 0 ? mag(num - exact)/mag(exact) : 0;
                const scalar modeErr =
                    exact != 0 ? residual/(mag(exact)*Foam::sqrt(p2)) : 0;
                const bool referenceMode = mag(theta - referenceTheta) <= SMALL;

                // Errors relative to the reference magnitude always apply;
                // relative to the mode's own eigenvalue only where it is
                // nonzero
                const bool pass =
                    std::isfinite(num) && std::isfinite(modeErrRef)
                 && rayleighErr <= tol && modeErrRef <= tol
                 && (exact == 0 || (relative <= tol && modeErr <= tol));

                failures += !pass;

                // Empty CSV fields where a value does not apply
                results << normalise << ',';
                if (normalise)
                {
                    results << rStar;
                }
                results
                    << ',' << names[modeli] << ',' << type << ','
                    << spectralRelation(type, modelDict, normalise) << ',';

                const label order = nLaplacianFactors(type, modelDict);
                if (order > 0)
                {
                    results << order;
                }
                results << ',' << scale << ',';
                if (power >= 0)
                {
                    results << power;
                }
                results
                    << ',' << csvQuote(parameters(modelDict)).c_str() << ','
                    << grid.N[0] << ',' << grid.N[1] << ',' << grid.N[2] << ','
                    << grid.h.x() << ',' << grid.h.y() << ',' << grid.h.z()
                    << ',' << modei << ','
                    << k[0] << ',' << k[1] << ',' << k[2] << ','
                    << theta.x()/pi << ',' << theta.y()/pi << ','
                    << theta.z()/pi << ','
                    << exact << ',' << num << ',' << ref << ',';
                if (exact != 0)
                {
                    results << relative;
                }
                results << ',' << rayleighErr << ',';
                if (exact != 0)
                {
                    results << modeErr;
                }
                results
                    << ',' << modeErrRef << ',' << referenceMode << ','
                    << pass << '\n';

                Info<< names[modeli]
                    << " (" << k[0] << ',' << k[1] << ',' << k[2] << ") "
                    << theta/pi << ' ' << exact << ' ' << num << ' ';
                if (exact != 0)
                {
                    Info<< relative << ' ' << modeErr;
                }
                else
                {
                    Info<< "n/a n/a (null residual=" << modeErrRef << ')';
                }
                Info<< ' ' << (pass ? "PASS" : "FAIL") << nl;

                if (optionFound(args, "writeFields"))
                {
                    const word suffix
                    (
                        word(normalise ? "normalised_" : "raw_")
                      + names[modeli] + "_mode" + Foam::name(modei)
                    );

                    volScalarField outputP
                    (
                        IOobject
                        (
                            "fourierP_" + suffix,
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        p
                    );
                    volScalarField outputS
                    (
                        IOobject
                        (
                            "fourierS_" + suffix,
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        stab
                    );
                    outputP.write();
                    outputS.write();
                }

                // Secondary check: the response is linear in a constant face
                // weight gamma, and a null weight gives exactly the response
                // with gamma = 1
                if (gammaCheck > 0)
                {
                    const scalarField scaled
                    (
                        primitiveField(model->cellScalar(&gammaC, true))
                    );
                    const scalar linearErr =
                        gMax(mag(scaled - gammaCheck*S))/ref;
                    const scalarField unweighted
                    (
                        primitiveField(model->cellScalar(nullptr, true))
                    );
                    const bool equal = std::memcmp
                    (
                        &unweighted[0],
                        &S[0],
                        S.size()*sizeof(scalar)
                    ) == 0;
                    const bool gammaPass =
                        std::isfinite(linearErr)
                     && linearErr <= gammaTol
                     && equal;

                    failures += !gammaPass;

                    gamma << normalise << ',';
                    if (normalise)
                    {
                        gamma << rStar;
                    }
                    gamma
                        << ',' << names[modeli] << ',' << modei << ','
                        << gammaCheck << ','
                        << linearErr << ',' << equal << ',' << gammaPass
                        << '\n';

                    if (!gammaPass)
                    {
                        Info<< "FAIL gamma check " << names[modeli]
                            << " mode " << modei << nl;
                    }
                }
            }

            // Pairs of models with the same operator must agree to pairTol;
            // the pairs that share code must agree bitwise
            forAll(pairs, pairi)
            {
                const label a = modelIndex(names, pairs[pairi].first());
                const label b = modelIndex(names, pairs[pairi].second());

                if (a < 0 || b < 0)
                {
                    FatalErrorInFunction
                        << "Unknown model in crossChecks " << pairs[pairi]
                        << exit(FatalError);
                }

                const scalar error =
                    gMax(mag(numerical[a] - numerical[b]))
                   /max(refs[a], refs[b]);
                const bool equal =
                    numerical[a].size() == numerical[b].size()
                 && std::memcmp
                    (
                        &numerical[a][0],
                        &numerical[b][0],
                        numerical[a].size()*sizeof(scalar)
                    ) == 0;
                const bool requireBitwise = sharedImplementation
                (
                    models.subDict(names[a]),
                    models.subDict(names[b])
                );
                const bool pass =
                    std::isfinite(error) && error <= pairTol
                 && (!requireBitwise || equal);

                failures += !pass;

                cross << normalise << ',';
                if (normalise)
                {
                    cross << rStar;
                }
                cross
                    << ',' << names[a] << ',' << names[b] << ',' << modei
                    << ',' << theta.x()/pi << ',' << theta.y()/pi << ','
                    << error << ',' << equal << ',' << pass << '\n';

                Info<< "pair " << names[a] << '/' << names[b]
                    << " mode " << modei
                    << " maxAbsDiff/ref=" << error << ' '
                    << (pass ? "PASS" : "FAIL") << nl;
            }
        }
    }

    if (optionFound(args, "symbolSweep"))
    {
        const bool written = writeSymbolExport
        (
            outDir,
            resolution,
            runTime,
            grid,
            models,
            names,
            tol,
            pairTol,
            gammaTol
        );

        if (!written)
        {
            ++failures;
        }
    }

    const bool resultsWritten = closeCsv(results);
    const bool crossWritten = closeCsv(cross);
    const bool gammaWritten = closeCsv(gamma);

    if (!resultsWritten || !crossWritten || !gammaWritten)
    {
        Info<< "FAIL: could not write verification CSVs" << nl;
        ++failures;
    }

    Info<< "Results: " << outDir << nl
        << "Failures: " << failures << nl
        << "End" << endl;

    return failures ? 1 : 0;
}


// ************************************************************************* //
