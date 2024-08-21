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

\*----------------------------------------------------------------------------*/

#include "sphericalCavityAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphericalCavityAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sphericalCavityAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::sphericalCavityAnalyticalSolution::KroneckerDelta
(
    const int i, const int j
) const
{
    return (i == j) ? 1.0 : 0.0;
}

Foam::scalar Foam::sphericalCavityAnalyticalSolution::calculateDispComponent
(
    const int i,
    const scalar nu,
    const scalar sigma_0,
    const scalar E,
    const scalar a,
    const vector& x
) const
{
    const scalar R = mag(x);
    const scalar x_3 = x[vector::Z];

    const scalar term1 =
        2
      + (5*(5 - 4*nu)*pow(a, 3))/((7 - 5*nu)*pow(R, 3))
      + (6*pow(a, 5))/((7 - 5*nu)*pow(R, 5));
    const scalar term2 =
        (-2*nu)/(1 + nu)
      + (5*nu - 6)*pow(a, 3)/((7 - 5*nu)*pow(R, 3))
      + 3*pow(a, 5)/((7 - 5*nu)*pow(R, 5));
    const scalar common_factor = (1 + nu)*sigma_0/(2*E);

    return
        common_factor
       *(
           term1*x_3*KroneckerDelta(i, 2)
         + term2*x[i]*(1 - 5*pow(x_3, 2)/pow(R, 2))
       );
}


Foam::symmTensor
Foam::sphericalCavityAnalyticalSolution::calculateSigma
(
    const scalar T,
    const scalar nu,
    const scalar a,
    const vector& x
) const
{
    // Stress formulae from Southwell et al, 1926, On the concentration of
    // stress in the neighbourhood of a small spherical flaw; and on the
    // propagation of fatigue fractures in Statistically Isotropic materials

    // Polar radial coordinate
    // Take care: in general, R != r
    const scalar R = mag(x);

    // Cylindrical radial coordinate
    const scalar r = sqrt(sqr(x.x()) + sqr(x.y()));

    // Cylindrical axial coordinate
    const scalar z = x.z();

    // Initialise stress tensor
    // We will place rr in xx, thetaTheta in yy, zz in zz, and zr in xz
    // We will late rotate this stress tensor to the Cartesian coordinate system
    symmTensor sigma = symmTensor::zero;

    // Radial (rr) component
    sigma.xx() =
        (T/(14 - 10*nu))*(pow3(a)/pow3(R))
       *(
            9 - 15*nu - 12*(sqr(a)/sqr(R))
          - (sqr(r)/sqr(R))*(72 - 15*nu - 105*(sqr(a)/sqr(R)))
          + 15*(pow4(r)/pow4(R))*(5 - 7*(sqr(a)/sqr(R)))
       );

    // Hoop (thetaTheta) component
    sigma.yy() =
        (T/(14 - 10*nu))*(pow3(a)/pow3(R))
       *(
            9 - 15*nu - 12*(sqr(a)/sqr(R))
          - 15*(sqr(r)/sqr(R))*(1 - 2*nu - (sqr(a)/sqr(R)))
       );

    // Axial (zz) component
    sigma.zz() =
        T
       *(
            1
          - (1/(14 - 10*nu))*(pow3(a)/pow3(R))
           *(
                38 - 10*nu - 24*sqr(a)/sqr(R)
              - (sqr(r)/sqr(R))*(117 - 15*nu - 120*sqr(a)/sqr(R))
              + 15*(pow4(r)/pow4(R))*(5 - 7*sqr(a)/sqr(R))
            )
        );

    // Radial-axial shear (rz) component
    sigma.xz() =
        (T/(14 - 10*nu))*(pow3(a)*z*r/pow5(R))
       *(
          - 3*(19 - 5*nu) + 60*(sqr(a)/sqr(R))
          + 15*(sqr(r)/sqr(R))*(5 - 7*(sqr(a)/sqr(R)))
        );

    return sigma;
}


// Foam::scalar
// Foam::sphericalCavityAnalyticalSolution::calculate_sigma_ij_over_sigma_0
// (
//     const int i,
//     const int j,
//     const scalar nu,
//     const scalar a,
//     const vector& x
// ) const
// {
//     const scalar R = mag(x);
//     const scalar x_3 = x[vector::Z];

//     const scalar delta_ij = KroneckerDelta(i, j);
//     const scalar delta_i3 = KroneckerDelta(i, 2);
//     const scalar delta_j3 = KroneckerDelta(j, 2);

//     const scalar term1 =
//         (3*pow(a, 3)/(2*(7 - 5*nu)*pow(R, 3)))
//        *(3 - 5*nu - 4*pow(a, 2)/pow(R, 2))*delta_ij;
//     const scalar term2 =
//         (3*pow(a, 3)*x[i]*x[j]/(2*(7 - 5*nu)*pow(R, 5)))
//        *(6 - 5*nu - 5*pow(a, 2)/pow(R, 2) + 10*pow(x_3, 2)/pow(R, 2));
//     const scalar term3 =
//         (delta_i3*delta_j3/(7 - 5*nu))
//        *((7 - 5*nu) + 5*(1 - 2*nu)*pow(a, 3)/pow(R, 3) + 3*pow(a, 5)/pow(R, 5));
//     const scalar term4 =
//         (-15*pow(a, 3)*x_3*(x[j]*delta_i3 + x[i]*delta_j3)/((7 - 5*nu)*pow(R, 5)))
//        *(pow(a, 2)/pow(R, 2) - nu);

//     return term1 + term2 + term3 + term4;
// }


Foam::vector
Foam::sphericalCavityAnalyticalSolution::sphericalCavityDisplacement
(
    const vector& C, const symmTensor& sigma
)
{
    // Calculate displacement
    vector disp = vector::zero;
    for (int i = 0; i < 3; ++i)
    {
        disp[i] = calculateDispComponent(i, nu_, T0_, E_, cavityR_, C);
    }

    return disp;
}


Foam::symmTensor
Foam::sphericalCavityAnalyticalSolution::sphericalCavityStress
(
    const vector& C
)
{
    // Calculate the stress in cylindrical coordinates
    const symmTensor sigmaCylindrical = calculateSigma(T0_, nu_, cavityR_, C);

    // Rotate the cylindrical stress to Cartesian coordinates 
    const coordinateSystem cs("cylindrical", C, vector(0, 0, 1), C/mag(C));
    const symmTensor sigmaCartesian = transform(cs.R(), sigmaCylindrical);

    return sigmaCartesian;
}


bool Foam::sphericalCavityAnalyticalSolution::writeData()
{
    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Lookup the point mesh
    const pointMesh& pMesh = mesh.lookupObject<pointMesh>("pointMesh");

    // Cell-centres coordinates
    const volVectorField& C = mesh.C();
    const vectorField& CI = C;

    // Point coordinates
    const pointField& points = mesh.points();

    if (gMin(mag(points)) < -SMALL)
    {
        FatalErrorIn("bool Foam::sphericalCavityAnalyticalSolution::writeData()")
            << "The cavity should be centred on the origin!" << endl
            << "gMin(mag(points)) = " << gMin(mag(points))
            << abort(FatalError);
    }

    // Cell analytical fields
    if (cellDisplacement_ || cellStress_)
    {
        // Analytical stress field
        volSymmTensorField analyticalStress
        (
            IOobject
            (
                "analyticalCellStress",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
            "calculated"
        );

        // Analytical displacement field
        volVectorField analyticalD
        (
            IOobject
            (
                "analyticalD",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimLength, vector::zero),
            "calculated"
        );

        symmTensorField& sI = analyticalStress;
        vectorField& aDI = analyticalD;

        forAll(sI, cellI)
        {
            if (cellStress_)
            {
                sI[cellI] = sphericalCavityStress(CI[cellI]);
            }

            if (cellDisplacement_)
            {
                aDI[cellI] = sphericalCavityDisplacement(CI[cellI], sI[cellI]);
            }
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
                vectorField& aDP = analyticalD.boundaryFieldRef()[patchI];
#else
                symmTensorField& sP = analyticalStress.boundaryField()[patchI];
                vectorField& aDP = analyticalD.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                    if (cellStress_)
                    {
                        sP[faceI] = sphericalCavityStress(CP[faceI]);
                    }

                    if (cellDisplacement_)
                    {
                        aDP[faceI] =
                            sphericalCavityDisplacement(CP[faceI], sP[faceI]);
                    }
                }
            }
        }

        // Write out the cell analytical field
        if (cellStress_)
        {
            Info<< "Writing analyticalCellStress field"
                << nl << endl;
            analyticalStress.write();
        }

        if (cellDisplacement_)
        {
            Info<< "Writing analyticalD field"
                << nl << endl;
            analyticalD.write();
        }


        if (cellStress_ && mesh.foundObject<volSymmTensorField>("sigma"))
        {
            const volSymmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>("sigma");

            const volSymmTensorField diff
            (
                "cellStressDifference", analyticalStress - sigma
            );
            Info<< "Writing cellStressDifference field" << endl;
            diff.write();

            for (int cmpt = 0; cmpt < pTraits<symmTensor>::nComponents; cmpt++)
            {
                // Only calculate for XX, XY and ZZ
                if (cmpt != 0 && cmpt != 1 && cmpt != 3)
                {
                    continue;
                }

                const symmTensorField& diffI = diff;
                const scalarField diffIcmptI(diffI.component(cmpt));

                Info<< "    Component: " << cmpt << endl;
                Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                    << "    " << gAverage(mag(diffIcmptI))
                    << " " << Foam::sqrt(gAverage(magSqr(diffIcmptI)))
                    << " " << gMax(mag(diffIcmptI))
                    << nl << endl;
            }
        }

        if (cellDisplacement_ && mesh.foundObject<volVectorField>("D"))
        {
            const volVectorField& D =
                mesh.lookupObject<volVectorField>("D");

            const volVectorField diff
            (
                "DDifference", analyticalD - D
            );
            Info<< "Writing DDifference field" << endl;
            diff.write();

            const vectorField& diffI = diff;
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << nl << endl;
        }
    }

    // Point analytical fields
    if (pointDisplacement_ || pointStress_)
    {
        pointSymmTensorField analyticalStress
        (
            IOobject
            (
                "analyticalPointStress",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        );

        pointVectorField analyticalD
        (
            IOobject
            (
                "analyticalPointD",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimLength, vector::zero)
        );

        symmTensorField& sI = analyticalStress;
        vectorField& aDI = analyticalD;

        forAll(sI, pointI)
        {
            if (pointStress_)
            {
                sI[pointI] = sphericalCavityStress(points[pointI]);
            }

            if (pointDisplacement_)
            {
                aDI[pointI] =
                    sphericalCavityDisplacement(points[pointI], sI[pointI]);
            }
        }

        // Write point analytical fields
        if (pointStress_)
        {
            Info<< "Writing analyticalPointStress"
                << nl << endl;
            analyticalStress.write();
        }

        if (pointDisplacement_)
        {
            Info<< "Writing analyticalPointDisplacement"
                << nl << endl;
            analyticalD.write();
        }

        if
        (
            pointDisplacement_
         && mesh.foundObject<pointVectorField>("pointD")
        )
        {
            const pointVectorField& pointD =
                mesh.lookupObject<pointVectorField>("pointD");

            const pointVectorField diff
            (
                "pointDDifference", analyticalD - pointD
            );
            Info<< "Writing pointDDifference field" << endl;
            diff.write();

            const vectorField& diffI = diff;
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << nl << endl;
        }
    }

    return true;
}

//** * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphericalCavityAnalyticalSolution::sphericalCavityAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    cellDisplacement_
    (
        dict.lookupOrDefault<Switch>("cellDisplacement", true)
    ),
    pointDisplacement_
    (
        dict.lookupOrDefault<Switch>("pointDisplacement", true)
    ),
    cellStress_
    (
        dict.lookupOrDefault<Switch>("cellStress", true)
    ),
    pointStress_
    (
        dict.lookupOrDefault<Switch>("pointStress", true)
    ),
    T0_(readScalar(dict.lookup("farFieldTractionZ"))),
    cavityR_(readScalar(dict.lookup("cavityRadius"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu")))
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (cavityR_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "cavityRadius should be greater than 0!"
            << abort(FatalError);
    }

    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sphericalCavityAnalyticalSolution::start()
{
    return true;
}


#if FOAMEXTEND
    bool Foam::sphericalCavityAnalyticalSolution::execute(const bool forceWrite)
#else
    bool Foam::sphericalCavityAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::sphericalCavityAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::sphericalCavityAnalyticalSolution::write()
{
    return false;
}
#endif

// ************************************************************************* //
