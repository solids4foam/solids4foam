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

#include "fvCFD.H"
#include "stabilisationModel.H"
#include "compatibilityFunctions.H"
#include "cartesianMeshInfo.H"
#include "fourierSymbols.H"
#include "OStringStream.H"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

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

static label modelIndex(const wordList& names, const word& name)
{
    forAll(names, i)
    {
        if (names[i] == name) return i;
    }
    return -1;
}

static bool splitSpectralShape(const word& type)
{
    return type == "diffStencilLaplacian" || type == "RhieChow";
}

static label paperOrder(const word& type, const dictionary& modelDict)
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
     && (type == "laplacian"
      || (type == "generalisedEvenOrderLaplacian"
       && readLabel(modelDict.lookup("laplacianPower")) == 0))
    )
    {
        return "legacyM1";
    }
    return "coupledRepeatedLaplacian";
}

//- Record arbitrary model settings, including nested sub-dictionaries.
static std::string parameters(const dictionary& dict, const std::string prefix = "")
{
    std::string result;
    const wordList keys(dict.toc());
    forAll(keys, i)
    {
        const std::string key(prefix + keys[i].c_str());
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
                if (j) value << ' ';
                value << tokens[j];
            }
            result += key + "=" + value.str().c_str() + ";";
        }
    }
    return result;
}

static std::string csvQuote(const std::string& value)
{
    std::string result("\"");
    for (const char c : value)
    {
        if (c == '"') result += '"';
        result += c;
    }
    return result + "\"";
}

static void metadata
(
    std::ofstream& out,
    const Time& runTime,
    const scalar tolerance,
    const scalar pairTolerance,
    const scalar gammaTolerance
)
{
    out << std::setprecision(17)
        << "# OpenFOAM=" << (std::getenv("WM_PROJECT_VERSION") ? std::getenv("WM_PROJECT_VERSION") : "unknown") << "\n# case=" << runTime.path()
        << "\n# gamma=1 (rAUf dimensions area/pressure)\n# tolerance=" << tolerance
        << "\n# pairTolerance=" << pairTolerance
        << "\n# gammaTolerance=" << gammaTolerance << '\n';
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
    argList::addOption("dict", "name", "Dictionary in system (default stabilisationFourierDict)");
    argList::addBoolOption("writeFields", "Write Fourier pressure and production stabilisation fields");
    argList::addBoolOption("symbolSweep", "Export analytical cuts and 2-D maps");
    argList::addOption("gammaLinearityCheck", "scalar", "Override secondary constant gamma (0 disables)");
    argList::addBoolOption("constructOnly", "Construct configured models, then exit");
#endif

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    word dictName("stabilisationFourierDict");
    readOption(args, "dict", dictName);
    IOdictionary dict
    (
        IOobject(dictName, runTime.system(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    const scalar tol = dict.lookupOrDefault<scalar>("tolerance", 1e-9);
    const scalar pairTol = dict.lookupOrDefault<scalar>("pairTolerance", 1e-12);
    const scalar gammaTol = dict.lookupOrDefault<scalar>("gammaTolerance", 1e-14);
    const Switch normalisedVerification
    (
        dict.lookupOrDefault<Switch>("normalisedVerification", false)
    );
    const label rStar = normalisedVerification
      ? readLabel(dict.lookup("referenceNyquistDirections"))
      : -1;
    scalar gammaCheck = dict.lookupOrDefault<scalar>("gammaLinearityCheck", 0.37);
    readOption(args, "gammaLinearityCheck", gammaCheck);
    if
    (
        !std::isfinite(tol) || tol <= 0
     || !std::isfinite(pairTol) || pairTol <= 0
     || !std::isfinite(gammaTol) || gammaTol <= 0
     || !std::isfinite(gammaCheck) || gammaCheck < 0
    )
    {
        FatalErrorInFunction << "Invalid tolerance or constant gamma"
            << exit(FatalError);
    }
    const cartesianMeshInfo grid(mesh);
    const scalar pi = std::acos(-1.0);
    const dictionary& models = dict.subDict("models");
    const wordList names(models.toc());
    if (optionFound(args, "constructOnly"))
    {
        forAll(names, modeli)
        {
            const dictionary& modelDict(models.subDict(names[modeli]));
            autoPtr<stabilisationModel> model
            (
                stabilisationModel::New
                (
                    mesh,
                    modelDict,
                    dimPressure/dimLength
                )
            );
            Info<< "Constructed stabilisation model " << names[modeli] << nl;
        }
        Info<< "End" << endl;
        return 0;
    }

    if (optionFound(args, "writeFields") && dict.found("fieldOutput"))
    {
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
            outputDict.lookupOrDefault<word>
            (
                "stabilisationPrefix",
                "fourierS"
            )
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
            const vector thetaByPi(outputModes.lookup(modeName));
            vector theta(vector::zero);
            label k[3];
            for (direction d = 0; d < 3; ++d)
            {
                const scalar discrete = grid.N[d]*thetaByPi[d]/2;
                if
                (
                    !std::isfinite(discrete)
                 || discrete < 0
                 || discrete > grid.N[d]/2.0
                )
                {
                    FatalErrorInFunction
                        << "Field-output mode " << modeName
                        << " must lie in [0,pi]" << exit(FatalError);
                }
                k[d] = std::lround(discrete);
                if
                (
                    mag(discrete - k[d]) > 1e-12
                 || (mesh.geometricD()[d] == -1 && thetaByPi[d] != 0)
                )
                {
                    FatalErrorInFunction
                        << "Field-output mode " << modeName
                        << " is incompatible with mesh periodicity"
                        << exit(FatalError);
                }
                theta[d] = 2*pi*k[d]/grid.N[d];
            }

            forAll(primitiveField(p), celli)
            {
                scalar phase = 0;
                for (direction d = 0; d < 3; ++d)
                {
                    phase += theta[d]*grid.index[d][celli];
                }
                primitiveFieldRef(p)[celli] = std::cos(phase);
            }
            p.correctBoundaryConditions();

            const word pName(pressurePrefix + "_" + modeName);
            volScalarField pMode
            (
                IOobject
                (
                    pName,
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
                const dictionary& modelDict =
                    models.subDict(names[modeli]);

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
                const volScalarField& stab =
                    model->cellScalar(&rAUf, true);
                const scalarField& S = primitiveField(stab);
                const scalar num = gSum(mesh.V()*pressure*S)/p2;
                const scalar residual = std::sqrt
                (
                    gSum(mesh.V()*sqr(S - num*pressure))
                );
                const scalar modeErr = residual
                  / max(mag(num)*std::sqrt(p2), SMALL);
                failures += !std::isfinite(num) || modeErr > tol;

                const word sName
                (
                    stabilisationPrefix + "_"
                  + names[modeli] + "_" + modeName
                );
                volScalarField outputS
                (
                    IOobject
                    (
                        sName,
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
            << runTime.path()/runTime.timeName() << nl
            << "Failures: " << failures << nl << "End" << endl;
        return failures ? 1 : 0;
    }

    const List<vector> modes(dict.lookup("modes"));
    const List<Pair<word>> pairs(dict.lookup("crossChecks"));
    if (names.empty() || modes.empty())
    {
        FatalErrorInFunction << "Models and modes must not be empty"
            << exit(FatalError);
    }

    const word resolution
    (
        Foam::name(grid.N[0]) + "x" + Foam::name(grid.N[1]) + "x" + Foam::name(grid.N[2])
    );
    const fileName outDir(runTime.path()/"postProcessing"/"stabilisationFourierCheck");
    mkDir(outDir);
    std::ofstream results((outDir/("results_" + resolution + ".csv")).c_str());
    std::ofstream cross((outDir/("crossChecks_" + resolution + ".csv")).c_str());
    std::ofstream gamma((outDir/("gammaChecks_" + resolution + ".csv")).c_str());
    metadata(results, runTime, tol, pairTol, gammaTol);
    metadata(cross, runTime, tol, pairTol, gammaTol);
    metadata(gamma, runTime, tol, pairTol, gammaTol);
    results << "normalise,referenceNyquistDirections,model,type,spectralRelation,paperOrder,"
        << "scaleFactor,laplacianPower,params,Nx,Ny,Nz,hx,hy,hz,modeIndex,kx,ky,kz,"
        << "thetaX_over_pi,thetaY_over_pi,thetaZ_over_pi,lambda_exact,lambda_num,lambda_ref,"
        << "relativeEigenvalueError,rayleighErr,eigenmodeErr,eigenmodeErrRef,referenceMode,pass\n";
    cross << "normalise,referenceNyquistDirections,modelA,modelB,modeIndex,"
        << "thetaX_over_pi,thetaY_over_pi,maxAbsDiffOverRef,bitwiseEqual,pass\n";
    gamma << "normalise,referenceNyquistDirections,model,modeIndex,gamma,"
        << "maxAbsDiffOverRef,nullGammaBitwiseEqual,pass\n";

    wordList patchTypes(mesh.boundary().size());
    forAll(patchTypes, patchi) patchTypes[patchi] = mesh.boundary()[patchi].type();
    volScalarField p
    (
        IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("zero", dimPressure, 0), patchTypes
    );
    surfaceScalarField rAUf
    (
        IOobject("rAUf", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("one", dimArea/dimPressure, 1)
    );
    surfaceScalarField gammaC
    (
        IOobject("gammaC", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("c", dimArea/dimPressure, gammaCheck)
    );
    const scalarField& V = mesh.V();
    label failures = 0;
    const label nCampaigns = normalisedVerification ? 2 : 1;
    for (label campaigni = 0; campaigni < nCampaigns; ++campaigni)
    {
        const bool normalise = campaigni == 1;
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
        label k[3];
        vector theta(vector::zero);
        for (direction d = 0; d < 3; ++d)
        {
            const scalar discrete = grid.N[d]*modes[modei][d]/2;
            if (!std::isfinite(discrete) || discrete < 0 || discrete > grid.N[d]/2.0)
            {
                FatalErrorInFunction << "Mode must lie in [0,pi]"
                    << exit(FatalError);
            }
            k[d] = std::lround(discrete);
            if
            (
                mag(discrete - k[d]) > 1e-12
             || (mesh.geometricD()[d] == -1 && modes[modei][d] != 0)
            )
            {
                FatalErrorInFunction << "Mode " << modes[modei]
                    << " is incompatible with mesh periodicity"
                    << exit(FatalError);
            }
            theta[d] = 2*pi*k[d]/grid.N[d];
        }
        forAll(primitiveField(p), celli)
        {
            scalar phase = 0;
            for (direction d = 0; d < 3; ++d)
            {
                phase += theta[d]*grid.index[d][celli];
            }
            // Index origin avoids the identically-zero Nyquist cosine that
            // would result from evaluating cos(pi*(i+1/2)).
            primitiveFieldRef(p)[celli] = std::cos(phase);
        }
        p.correctBoundaryConditions();
        volVectorField gradP(fvc::grad(p));
        gradP.correctBoundaryConditions();
        const scalarField pressure(primitiveField(p));
        const scalar p2 = gSum(V*sqr(pressure));
        if (p2 <= SMALL)
        {
            FatalErrorInFunction << "Degenerate Fourier mode" << exit(FatalError);
        }
        List<scalarField> numerical(names.size());
        scalarField refs(names.size());
        forAll(names, modeli)
        {
            dictionary modelDict(models.subDict(names[modeli]));
            if (normalise)
            {
                modelDict.set("normalise", Switch(true));
                modelDict.set("referenceNyquistDirections", rStar);
            }
            const word type(modelDict.lookup("type"));
            const scalar scale = type == "laplacian"
              ? modelDict.lookupOrDefault<scalar>("scaleFactor", 1)
              : readScalar(modelDict.lookup("scaleFactor"));
            const label power = type == "generalisedEvenOrderLaplacian"
              ? readLabel(modelDict.lookup("laplacianPower")) : -1;
            const scalar exact = expectedEigenvalue(modelDict, theta, grid.h, grid.N);
            const scalar ref = mag(expectedEigenvalue(modelDict, referenceTheta, grid.h, grid.N));
            if (!std::isfinite(ref) || ref <= 0)
            {
                FatalErrorInFunction << "Verification needs a finite nonzero reference response"
                    << exit(FatalError);
            }
            refs[modeli] = ref;

            // Numerical result: exclusively the production runtime-selected
            // model, its updateScalar and its intensive cellScalar interface.
            autoPtr<stabilisationModel> model
            (
                stabilisationModel::New(mesh, modelDict, dimPressure/dimLength)
            );
            model->updateScalar(p, &gradP);
            const volScalarField& stab = model->cellScalar(&rAUf, true);
            numerical[modeli] = primitiveField(stab);
            const scalarField& S = numerical[modeli];
            const scalar num = gSum(V*pressure*S)/p2;
            const scalar residual = std::sqrt(gSum(V*sqr(S - exact*pressure)));
            const scalar rayleighErr = mag(num - exact)/ref;
            const scalar modeErrRef = residual/(ref*std::sqrt(p2));
            const scalar relative = exact != 0 ? mag(num - exact)/mag(exact) : 0;
            const scalar modeErr = exact != 0 ? residual/(mag(exact)*std::sqrt(p2)) : 0;
            const bool referenceMode = mag(theta - referenceTheta) <= SMALL;
            const bool pass = std::isfinite(num) && std::isfinite(modeErrRef)
                && rayleighErr <= tol && modeErrRef <= tol
                && (exact == 0 || (relative <= tol && modeErr <= tol));
            failures += !pass;
            results << normalise << ',';
            if (normalise) results << rStar;
            results << ',' << names[modeli] << ',' << type << ','
                << spectralRelation(type, modelDict, normalise) << ',';
            const label order = paperOrder(type, modelDict);
            if (order > 0) results << order;
            results << ',' << scale << ',';
            if (power >= 0) results << power;
            results << ',' << csvQuote(parameters(modelDict)) << ','
                << grid.N[0] << ',' << grid.N[1] << ',' << grid.N[2] << ','
                << grid.h.x() << ',' << grid.h.y() << ',' << grid.h.z() << ','
                << modei << ',' << k[0] << ',' << k[1] << ',' << k[2] << ','
                << theta.x()/pi << ',' << theta.y()/pi << ',' << theta.z()/pi << ','
                << exact << ',' << num << ',' << ref << ',';
            if (exact != 0) results << relative;
            results << ',' << rayleighErr << ',';
            if (exact != 0) results << modeErr;
            results << ',' << modeErrRef << ',' << referenceMode << ','
                << pass << '\n';
            Info<< names[modeli] << " (" << k[0] << ',' << k[1] << ',' << k[2] << ") "
                << theta/pi << ' ' << exact << ' ' << num << ' ';
            if (exact != 0) Info<< relative << ' ' << modeErr;
            else Info<< "n/a n/a (null residual=" << modeErrRef << ')';
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
                    IOobject("fourierP_" + suffix, runTime.timeName(), mesh,
                        IOobject::NO_READ, IOobject::NO_WRITE, false), p
                );
                volScalarField outputS
                (
                    IOobject("fourierS_" + suffix, runTime.timeName(), mesh,
                        IOobject::NO_READ, IOobject::NO_WRITE, false), stab
                );
                outputP.write();
                outputS.write();
            }
            if (gammaCheck > 0)
            {
                const scalarField scaled(primitiveField(model->cellScalar(&gammaC, true)));
                const scalar linearErr = gMax(mag(scaled - gammaCheck*S))/ref;
                const scalarField unweighted(primitiveField(model->cellScalar(nullptr, true)));
                const bool equal = std::memcmp
                (
                    &unweighted[0], &S[0], S.size()*sizeof(scalar)
                ) == 0;
                const bool gammaPass = std::isfinite(linearErr) && linearErr <= gammaTol && equal;
                failures += !gammaPass;
                gamma << normalise << ',';
                if (normalise) gamma << rStar;
                gamma << ',' << names[modeli] << ',' << modei << ','
                    << gammaCheck << ','
                    << linearErr << ',' << equal << ',' << gammaPass << '\n';
                if (!gammaPass) Info<< "FAIL gamma check " << names[modeli] << " mode " << modei << nl;
            }
        }
        forAll(pairs, pairi)
        {
            const label a = modelIndex(names, pairs[pairi].first());
            const label b = modelIndex(names, pairs[pairi].second());
            if (a < 0 || b < 0)
            {
                FatalErrorInFunction << "Unknown model in crossChecks " << pairs[pairi]
                    << exit(FatalError);
            }
            const scalar error = gMax(mag(numerical[a] - numerical[b]))/max(refs[a], refs[b]);
            const bool equal = numerical[a].size() == numerical[b].size()
             && std::memcmp
                (
                    &numerical[a][0],
                    &numerical[b][0],
                    numerical[a].size()*sizeof(scalar)
                ) == 0;
            const bool requireBitwise =
                (
                    names[a] == "laplacian" && names[b] == "gen0"
                )
             || (
                    names[a] == "diffStencil" && names[b] == "RhieChow"
                );
            const bool pass = std::isfinite(error) && error <= pairTol
                && (!requireBitwise || equal);
            failures += !pass;
            cross << normalise << ',';
            if (normalise) cross << rStar;
            cross << ',' << names[a] << ',' << names[b] << ',' << modei << ','
                << theta.x()/pi << ',' << theta.y()/pi << ',' << error << ','
                << equal << ',' << pass << '\n';
            Info<< "pair " << names[a] << '/' << names[b] << " mode " << modei
                << " maxAbsDiff/ref=" << error << ' ' << (pass ? "PASS" : "FAIL") << nl;
        }
    }
    }

    if (optionFound(args, "symbolSweep"))
    {
        std::ofstream sweep((outDir/("symbols_" + resolution + ".csv")).c_str());
        std::ofstream maps((outDir/("symbolMaps_" + resolution + ".csv")).c_str());
        metadata(sweep, runTime, tol, pairTol, gammaTol);
        metadata(maps, runTime, tol, pairTol, gammaTol);
        sweep << "model,cut,theta_over_pi,lambda_exact,hx,hy,scaleFactor\n";
        maps << "model,thetaX_over_pi,thetaY_over_pi,lambda_exact,hx,hy,scaleFactor\n";
        forAll(names, modeli)
        {
            const dictionary& md = models.subDict(names[modeli]);
            const scalar scale = md.lookupOrDefault<scalar>("scaleFactor", 1);
            for (label cut = 0; cut < 3; ++cut)
            {
                for (label i = 0; i <= 200; ++i)
                {
                    const scalar t = pi*i/200;
                    const vector theta(cut == 2 ? pi : t, cut == 0 ? 0 : t, 0);
                    sweep << names[modeli] << ',' << cut << ',' << t/pi << ','
                        << expectedEigenvalue(md, theta, grid.h, grid.N) << ','
                        << grid.h.x() << ',' << grid.h.y() << ',' << scale << '\n';
                }
            }
            for (label j = 0; j <= 64; ++j)
            {
                for (label i = 0; i <= 64; ++i)
                {
                    const vector theta(pi*i/64, pi*j/64, 0);
                    maps << names[modeli] << ',' << i/64.0 << ',' << j/64.0 << ','
                        << expectedEigenvalue(md, theta, grid.h, grid.N) << ','
                        << grid.h.x() << ',' << grid.h.y() << ',' << scale << '\n';
                }
            }
        }
        sweep.close();
        maps.close();
        if (!sweep || !maps) ++failures;
    }
    results.close();
    cross.close();
    gamma.close();
    if (!results || !cross || !gamma)
    {
        Info<< "FAIL: could not write verification CSVs" << nl;
        ++failures;
    }
    Info<< "Results: " << outDir << nl << "Failures: " << failures << nl << "End" << endl;
    return failures ? 1 : 0;
}

// ************************************************************************* //
