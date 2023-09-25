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
    const word modelType(dict.lookup("type"));

    Info<< "Selecting mechanical law " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = linGeomMechLawConstructorTable(modelType);

    if (!ctorPtr)
    {
        // Check if the user inadvertently specified a nonlinear law
        auto* ctorPtrOther = nonLinGeomMechLawConstructorTable(modelType);

        if (ctorPtrOther)
        {
            FatalIOErrorInFunction(dict)
                << "The mechanicalLaw " << modelType
                << " can only be used with a linear geometry solid model"
                << endl << endl
                << "Valid linearGeometry mechanicalLaws are : " << endl
                << linGeomMechLawConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorInLookup
            (
                dict,
                "mechanicalLaw",
                modelType,
                *linGeomMechLawConstructorTablePtr_
            )<< exit(FatalIOError);
        }
    }

#else
    linGeomMechLawConstructorTable::iterator cstrIter =
        linGeomMechLawConstructorTablePtr_->find(modelType);

    if (cstrIter == linGeomMechLawConstructorTablePtr_->end())
    {
        // Check if the user inadvertently specified a nonlinear law
        nonLinGeomMechLawConstructorTable::iterator nlgLawIter =
            nonLinGeomMechLawConstructorTablePtr_->find(modelType);

        if (nlgLawIter != nonLinGeomMechLawConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "mechanicalLaw::New(\n"
                "    const word& name,\n"
                "    const fvMesh& mesh,\n"
                "    dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "The mechanicalLaw " << modelType
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
                "    const fvMesh& mesh,\n"
                "    dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "Unknown mechanicalLaw type "
                << modelType << endl << endl
                << "Valid linearGeometry mechanicalLaws are : " << endl
                << linGeomMechLawConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<mechanicalLaw>(ctorPtr(name, mesh, dict, nonLinGeom));
}


autoPtr<mechanicalLaw> mechanicalLaw::NewNonLinGeomMechLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
{
    const word modelType(dict.lookup("type"));

    Info<< "Selecting mechanical law " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = nonLinGeomMechLawConstructorTable(modelType);

    if (!ctorPtr)
    {
        // Check if the user inadvertently specified a nonlinear law
        auto* ctorPtrOther = linGeomMechLawConstructorTable(modelType);

        if (ctorPtrOther)
        {
            FatalIOErrorInFunction(dict)
                << "The mechanicalLaw " << modelType
                << " can only be used with a non-linear geometry solid model"
                << endl << endl
                << "Valid nonLinearGeometry mechanicalLaws are : " << endl
                << nonLinGeomMechLawConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorInLookup
            (
                dict,
                "mechanicalLaw",
                modelType,
                *nonLinGeomMechLawConstructorTablePtr_
            ) << exit(FatalIOError);
        }
    }

    return autoPtr<mechanicalLaw>(ctorPtr(name, mesh, dict, nonLinGeom));

    #else
    nonLinGeomMechLawConstructorTable::iterator cstrIter =
        nonLinGeomMechLawConstructorTablePtr_->find(modelType);

    if (cstrIter == nonLinGeomMechLawConstructorTablePtr_->end())
    {
        // Check if the user inadvertently specified a nonlinear law
        linGeomMechLawConstructorTable::iterator lgLawIter =
            linGeomMechLawConstructorTablePtr_->find(modelType);

        if (lgLawIter != linGeomMechLawConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "mechanicalLaw::New(\n"
                "    const word& name,\n"
                "    const fvMesh& mesh,\n"
                "    dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "The mechanicalLaw " << modelType
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
                "    const fvMesh& mesh,\n"
                "    dictionary& dict,\n"
                "    const nonLinearGeometry::nonLinearType& nonLinGeom\n"
                ")",
                dict
            )   << "Unknown mechanicalLaw type "
                << modelType << endl << endl
                << "Valid nonLinearGeometry mechanicalLaws are : " << endl
                << nonLinGeomMechLawConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
    }

    auto* ctorPtr = cstrIter();
    #endif

    return autoPtr<mechanicalLaw>(ctorPtr(name, mesh, dict, nonLinGeom));

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
