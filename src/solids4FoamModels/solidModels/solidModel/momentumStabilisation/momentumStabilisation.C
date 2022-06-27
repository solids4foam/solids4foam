/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "momentumStabilisation.H"
#include "fvc.H"

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
#ifdef OPENFOAMESIORFOUNDATION
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

#ifdef OPENFOAMESIORFOUNDATION
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

#ifdef OPENFOAMESIORFOUNDATION
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
    }

    // Calculate stabilisation term
    if (method == "RhieChow")
    {
        result =
        (
           fvc::laplacian(gammaf, vf, "laplacian(DD,D)")
         - fvc::div(gammaf*mesh.Sf() & fvc::interpolate(gradVf))
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
    else if (method == "Laplacian")
    {
        result = fvc::laplacian(gammaf, vf);
    }
    else if (method != "none")
    {
        FatalErrorIn(type() + "::stabilisation() const")
            << "Unknown method = " << method << nl
            << "Methods are: none, RhieChow, JamesonSchmidtTurkel and Laplacian"
            <<  abort(FatalError);
    }

    return tresult;
}



// ************************************************************************* //
