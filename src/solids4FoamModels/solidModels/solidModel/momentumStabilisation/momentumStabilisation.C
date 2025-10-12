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

#include "momentumStabilisation.H"
#include "fvc.H"
#include "deltaVectors.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentumStabilisation, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumStabilisation::momentumStabilisation
(
    const dictionary& dict
)
:
    dict_(dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::momentumStabilisation::stabilisation
(
    const volVectorField& vf,
    const volTensorField& gradVf,
    const volScalarField& gamma
) const
{
    const fvMesh& mesh = vf.mesh();

    tmp<volVectorField> tresult
    (
        new volVectorField
        (
            IOobject
            (
                word(type() + "Stabilisation"),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zero",
                gamma.dimensions()*gradVf.dimensions()/dimLength,
                vector::zero
            )
        )
    );
#ifdef OPENFOAM_NOT_EXTEND
    volVectorField& result = tresult.ref();
#else
    volVectorField& result = tresult();
#endif

    // Lookup method
    const word method = word(dict_.lookup("type"));

    // Merge scale factor into gamma field
    volScalarField gammaMod("gammaStabilisation", gamma);
    scalarField& gammaModI = gammaMod;

    // Optional: cell zone scale factors
    if (method != "none")
    {
        if (dict_.lookupOrDefault<Switch>("useCellZones", false))
        {
            forAll(vf.mesh().cellZones(), czI)
            {
                const cellZone& cz = vf.mesh().cellZones()[czI];

                // Lookup the scale factor for this cell zone
                const scalar czScaleFac =
                    readScalar(dict_.lookup(cz.name() + "ScaleFactor"));

                forAll(cz, cI)
                {
                    const label cellI = cz[cI];
                    gammaModI[cellI] *= czScaleFac;
                }
            }

            // Correct processor boundaries
            gammaMod.correctBoundaryConditions();
        }
        else
        {
            // Read scale factor
            const scalar scaleFactor = readScalar(dict_.lookup("scaleFactor"));

            gammaMod *= scaleFactor;
        }
    }

    // Interpolate gamma to faces
    surfaceScalarField gammaf(fvc::interpolate(gammaMod, "interpolate(impK)"));

    // Specify different scale factor at material interfaces
    const scalar interfaceScaleFactor
    (
        dict_.lookupOrDefault<scalar>("interfaceScaleFactor", 0.01)
    );

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& gammaI = gamma;

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& gammafI = gammaf.ref();
#else
    scalarField& gammafI = gammaf.internalField();
#endif

    forAll(gammafI, faceI)
    {
        const scalar gOwn = gammaI[own[faceI]];
        const scalar gNei = gammaI[nei[faceI]];

        if (mag(gOwn - gNei) > SMALL)
        {
            // Harmonic
            // gammafI[faceI] = interfaceScaleFactor*(gOwn*gNei)/(gOwn + gNei);

            // Arithmetric average
            gammafI[faceI] = interfaceScaleFactor*0.5*(gOwn + gNei);
        }
    }

    // Correct processor patches
    forAll(gammaf.boundaryField(), patchI)
    {
        if (vf.boundaryField()[patchI].coupled())
        {
            const scalarField pif
            (
                gamma.boundaryField()[patchI].patchInternalField()
            );
            const scalarField pnf
            (
                gamma.boundaryField()[patchI].patchNeighbourField()
            );

#ifdef OPENFOAM_NOT_EXTEND
            scalarField& gammafP = gammaf.boundaryFieldRef()[patchI];
#else
            scalarField& gammafP = gammaf.boundaryField()[patchI];
#endif

            if (vf.boundaryField()[patchI].type() == "processor")
            {
                // Some faces may be on a material interface but some may
                // not
                forAll(gammafP, faceI)
                {
                    const scalar gOwn = pif[faceI];
                    const scalar gNei = pnf[faceI];

                    if (mag(gOwn - gNei) > SMALL)
                    {
                        // Harmonic
                        // gammafP[faceI] =
                        // interfaceScaleFactor*(gOwn*gNei)/(gOwn + gNei);

                        // Arithmetric average
                        gammafP[faceI] =
                            interfaceScaleFactor*0.5*(gOwn + gNei);
                    }
                }
            }
            else
            {
                gammafP = interfaceScaleFactor*0.5*(pif + pnf);
            }
        }
        else
        {
            // Set stabilisation to zero on non-coupled boundaries
#ifdef OPENFOAM_NOT_EXTEND
            gammaf.boundaryFieldRef()[patchI] = 0.0;
#else
            gammaf.boundaryField()[patchI] = 0.0;
#endif
        }
    }

    // Calculate stabilisation term
    if (method == "RhieChow")
    {
        result =
        (
           fvc::laplacian(gammaf, vf, "laplacian(DD,D)")
         - fvc::div(gammaf*(mesh.Sf() & fvc::interpolate(gradVf)))
        );
    }
    else if (method == "JamesonSchmidtTurkel")
    {
        result = -fvc::laplacian
        (
            vf.mesh().magSf(),
            fvc::laplacian(gammaf, vf, "JSTinner"),
            "JSTouter"
        );
    }
    else if (method == "alpha")
    {
        // This alpha scheme is adapted from: H Nishikawa, Y Nakashima, N
        // Watanabe, Effects of high-frequency damping on iterative convergence
        // of implicit viscous solver, Journal of Computational Physics, 2017,
        // 10.1016/j.jcp.2017.07.021.

        // Note that the scale factor `alpha` has already been rolled into the
        // gammaf field above

        const surfaceVectorField n(mesh.Sf()/mesh.magSf());
        const vectorField& nI = n.internalField();
        const vectorField deltaI(deltaVectors(mesh));
        const scalarField& magSfI = mesh.magSf().internalField();
        const vectorField& vfI = vf.internalField();
        const tensorField& gradVfI = gradVf.internalField();

        vectorField& resultI = result;

        forAll(gammafI, faceI)
        {
            const label ownCellID = own[faceI];
            const label neiCellID = nei[faceI];

            const vector& ownVf = vfI[ownCellID];
            const vector& neiVf = vfI[neiCellID];

            const tensor& ownGradVf = gradVfI[ownCellID];
            const tensor& neiGradVf = gradVfI[neiCellID];

            const vector& n = nI[faceI];
            const vector& d = deltaI[faceI];

            const vector extrapOwnVf = ownVf + 0.5*(d & ownGradVf);
            const vector extrapNeiVf = neiVf - 0.5*(d & neiGradVf);

            const scalar curMagSf = magSfI[faceI];
            const scalar curGamma = gammafI[faceI];

            // Note: gamma already contains the scale factor alpha
            const vector faceDamping =
                curMagSf*curGamma*(extrapNeiVf - extrapOwnVf)/mag(n & d);

            // The alpha scheme is essentially equivalent to Rhie-Chow; the
            // differences may lie in how the weights are calculated, which may
            // have consequences for distorted meshes: the gradVf terms are
            // scaled by mag(n & d) in the alpha scheme but not in the Rhie-Chow
            // scheme; for a highly non-orthogonal face, mag(n & d) will be
            // small, which increases the entire damping term for the alpha
            // scheme; whereas only the compact term is increased in the
            // Rhie-Chow approach.

            // Rhie-Chow (mid-point interpolation of gradients)
            // const vector faceDamping =
            //     curMagSf*curGamma
            //    *(
            //        (neiVf - ownVf)/(n & d)
            //      - 0.5*(n & (ownGradVf + neiGradVf))
            //    );

            resultI[ownCellID] += faceDamping;
            resultI[neiCellID] -= faceDamping;
        }

        if (Pstream::parRun())
        {
            notImplemented
            (
                "Parallel boundaries for alpha scheme have to be implemented"
            );
        }

        // Divide by the volume
        resultI /= mesh.V();
    }
    else if (method == "Laplacian")
    {
        result = fvc::laplacian(gammaf, vf);
    }
    else if (method != "none")
    {
        FatalErrorIn(type() + "::stabilisation() const")
            << "Unknown method = " << method << nl
            << "Methods are: none, RhieChow, JamesonSchmidtTurkel, alpha and "
            << "Laplacian" <<  abort(FatalError);
    }

    return tresult;
}



// ************************************************************************* //
