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

#include "mmsFsiSource.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "mmsFsiFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(mmsFsiSource, 0);
    addToRunTimeSelectionTable(option, mmsFsiSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::vectorField& Foam::fv::mmsFsiSource::source() const
{
    const scalar t = mesh_.time().value();
    const vectorField& C = mesh_.C().primitiveField();

    bool update = (t != sourceTime_) || (C.size() != sourceC_.size());

    if (!update)
    {
        forAll(C, cellI)
        {
            if (C[cellI] != sourceC_[cellI])
            {
                update = true;
                break;
            }
        }
    }

    if (update)
    {
        sourceTime_ = t;
        sourceC_ = C;
        source_.setSize(C.size());

        if (region_ == "fluid")
        {
            forAll(C, cellI)
            {
                source_[cellI] = mmsFsiFields::fluid(C[cellI], t).f;
            }
        }
        else
        {
            forAll(C, cellI)
            {
                source_[cellI] = mmsFsiFields::solid(C[cellI], t).b0;
            }
        }
    }

    return source_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::mmsFsiSource::mmsFsiSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    region_(coeffs_.get<word>("region")),
    sourceTime_(-GREAT),
    sourceC_(),
    source_()
{
    coeffs_.readEntry("fields", fieldNames_);

    if (region_ != "fluid" && region_ != "solid")
    {
        FatalIOErrorInFunction(dict)
            << "region must be fluid or solid, not " << region_
            << exit(FatalIOError);
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::mmsFsiSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    // The source has the dimensions of the equation per unit volume
    DimensionedField<vector, volMesh> S
    (
        IOobject
        (
            "mmsFsiSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<vector>(eqn.dimensions()/dimVolume, Zero)
    );
    S.field() = source();

    eqn += S;
}


void Foam::fv::mmsFsiSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup(eqn, fieldi);
}


bool Foam::fv::mmsFsiSource::read(const dictionary& dict)
{
    return fv::option::read(dict);
}


// ************************************************************************* //
