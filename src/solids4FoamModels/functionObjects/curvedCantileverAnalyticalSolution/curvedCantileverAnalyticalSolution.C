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

#include "curvedCantileverAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(curvedCantileverAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        curvedCantileverAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor Foam::curvedCantileverAnalyticalSolution::curvedCantileverStress
(
    const vector& C
)
{
    tensor sigma = tensor::zero;
    
    const scalar r = Foam::sqrt(Foam::sqr(C.x()) + Foam::sqr(C.y()));
    
    if (r < SMALL)
    {
        FatalErrorIn("Foam::symmTensor Foam::"
            "curvedCantileverAnalyticalSolution::curvedCantileverStress"
            "(const vector& C)")
            << "The beam radius in zero: this is not allowed! "
            << "The beam's center of curvature should be centered at the origin"
            << abort(FatalError);
    }
       
    const scalar theta = Foam::atan2(C.y(), C.x());
    
    const coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    const scalar& a = rInner_;
    const scalar& b = rOuter_;
    
    const scalar N = Foam::sqr(a) - Foam::sqr(b) 
        + (Foam::sqr(a) + Foam::sqr(b))*Foam::log(b/a);
        
    sigma.xx() =
        - (force_ / N) * Foam::sin(theta)
      * (
            r 
         + ((Foam::sqr(a)*Foam::sqr(b))/Foam::pow(r,3))
         - ((Foam::sqr(a)+Foam::sqr(b))/r)
        );
        
    sigma.yy() =
        - (force_ / N) * Foam::sin(theta)
      * (
            3*r 
         - ((Foam::sqr(a)*Foam::sqr(b))/Foam::pow(r,3))
         - ((Foam::sqr(a)+Foam::sqr(b))/r)
        );
        
    sigma.xy() =
        (force_ / N) * Foam::cos(theta)
      * (
            r 
         + ((Foam::sqr(a)*Foam::sqr(b))/Foam::pow(r,3))
         - ((Foam::sqr(a)+Foam::sqr(b))/r)
        );
        
    sigma.yx() = sigma.xy();
    
    // Transformation to Cartesian coordinate system
#ifdef OPENFOAMFOUNDATION
    sigma = ((cs.R().R() & sigma) & cs.R().R().T());
#else
    sigma = ((cs.R() & sigma) & cs.R().T());
#endif

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}

bool Foam::curvedCantileverAnalyticalSolution::writeData()
{
    Info<< "curvedCantileverAnalyticalSolution:" << nl
        << "Writing analytical solution fields"
        << endl;

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

    // Cell centre coordinates
    const volVectorField& C = mesh.C();
    const vectorField& CI = C;

    volSymmTensorField analyticalStress
    (
        IOobject
        (
           "analyticalStress",
           time_.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    symmTensorField& sI = analyticalStress;

    forAll(sI, cellI)
    {
        sI[cellI] = curvedCantileverStress(CI[cellI]);
    }

    forAll(analyticalStress.boundaryField(), patchI)
    {
        if (mesh.boundary()[patchI].type() != "empty")
        {
#ifdef OPENFOAMESIORFOUNDATION
            symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
#else
            symmTensorField& sP = analyticalStress.boundaryField()[patchI];
#endif
            const vectorField& CP = C.boundaryField()[patchI];

            forAll(sP, faceI)
            {
                sP[faceI] = curvedCantileverStress(CP[faceI]);
            }
        }
    }

    analyticalStress.write();

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvedCantileverAnalyticalSolution::curvedCantileverAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    rInner_(readScalar(dict.lookup("rInner"))),
    rOuter_(readScalar(dict.lookup("rOuter"))),
    force_(readScalar(dict.lookup("force"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu")))
{
   
    Info<< "Creating " << this->name() << " function object" << endl;

    if (rInner_ >= rOuter_)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "rInner should be less than rOuter!"
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

bool Foam::curvedCantileverAnalyticalSolution::start()
{
    return true;
}


#if FOAMEXTEND
bool Foam::curvedCantileverAnalyticalSolution::execute(const bool forceWrite)
#else
bool Foam::curvedCantileverAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::curvedCantileverAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::curvedCantileverAnalyticalSolution::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
