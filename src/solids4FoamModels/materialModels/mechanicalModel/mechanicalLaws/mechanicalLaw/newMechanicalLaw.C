/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

InClass
    mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<mechanicalLaw> mechanicalLaw::NewLinGeomMechLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
{
    const word mechTypeName(dict.lookup("type"));

    Info<< "Selecting mechanical law " << mechTypeName << endl;

    linGeomMechLawConstructorTable::iterator cstrIter =
        linGeomMechLawConstructorTablePtr_->find(mechTypeName);

    if (cstrIter == linGeomMechLawConstructorTablePtr_->end())
    {
        // Check if the user inadvertently specified a nonlinear law
        nonLinGeomMechLawConstructorTable::iterator nlgLawIter =
            nonLinGeomMechLawConstructorTablePtr_->find(mechTypeName);

        if (nlgLawIter != nonLinGeomMechLawConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "mechanicalLaw::New(\n"
                "    const word& name,\n"
                "    const fvMehs& mesh,\n"
                "    const dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "The mechanicalLaw " << mechTypeName
                << " can only be used with a non-linear geometry solid model"
                << endl << endl
                << "Valid linearGeometry mechanicalLaws are : " << endl
                << linGeomMechLawConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorIn
            (
                "mechanicalLaw::New(\n"
                "    const word& name,\n"
                "    const fvMehs& mesh,\n"
                "    const dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "Unknown mechanicalLaw type "
                << mechTypeName << endl << endl
                << "Valid linearGeometry mechanicalLaws are : " << endl
                << linGeomMechLawConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
    }

    return autoPtr<mechanicalLaw>(cstrIter()(name, mesh, dict, nonLinGeom));
}


autoPtr<mechanicalLaw> mechanicalLaw::NewNonLinGeomMechLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
{
    const word mechTypeName(dict.lookup("type"));

    Info<< "Selecting mechanical law " << mechTypeName << endl;

    nonLinGeomMechLawConstructorTable::iterator cstrIter =
        nonLinGeomMechLawConstructorTablePtr_->find(mechTypeName);

    if (cstrIter == nonLinGeomMechLawConstructorTablePtr_->end())
    {
        // Check if the user inadvertently specified a nonlinear law
        linGeomMechLawConstructorTable::iterator lgLawIter =
            linGeomMechLawConstructorTablePtr_->find(mechTypeName);

        if (lgLawIter != linGeomMechLawConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "mechanicalLaw::New(\n"
                "    const word& name,\n"
                "    const fvMehs& mesh,\n"
                "    const dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "The mechanicalLaw " << mechTypeName
                << " can only be used with a linear geometry solid model"
                << endl << endl
                << "Valid nonLinearGeometry mechanicalLaws are : " << endl
                << nonLinGeomMechLawConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorIn
            (
                "mechanicalLaw::New(\n"
                "    const word& name,\n"
                "    const fvMehs& mesh,\n"
                "    const dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "Unknown mechanicalLaw type "
                << mechTypeName << endl << endl
                << "Valid nonLinearGeometry mechanicalLaws are : " << endl
                << nonLinGeomMechLawConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
    }

    return autoPtr<mechanicalLaw>(cstrIter()(name, mesh, dict, nonLinGeom));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
