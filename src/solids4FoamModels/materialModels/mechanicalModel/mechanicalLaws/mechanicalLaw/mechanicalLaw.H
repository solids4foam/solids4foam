/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

Class
    Foam::mechanicalLaw

Description
    Mechanical constitutive law for the solidModels.

SourceFiles
    mechanicalLaw.C
    newMechanicalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalLaw_H
#define mechanicalLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class mechanicalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Mechanical law dictionary
        const dictionary dict_;

        //- Name of the baseMesh region
        word baseMeshRegionName_;

        //- Nonlinear geometry enumerator
        const nonLinearGeometry::nonLinearType nonLinGeom_;


    // Private Member Functions

        //- Disallow copy construct
        mechanicalLaw(const mechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalLaw&);


protected:

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return const reference to dict
        const dictionary& dict() const
        {
            return dict_;
        }

        //- If the case is plane stress
        bool planeStress() const;

        //- Return the base mesh region name
        const word& baseMeshRegionName() const
        {
            return baseMeshRegionName_;
        }

        //- Return const access to the base mesh
        const fvMesh& baseMesh() const
        {
            return mesh().db().lookupObject<fvMesh>
            (
                baseMeshRegionName()
            );
        }

        //- Lookup the enforceLinear Switch in the solidModel
        const Switch& enforceLinear() const;

        //- Does the solver use an incremental approach
        //  i.e. does it solve for DD as opposed to D
        bool incremental() const;

        //- Return nonlinear geometry enumerator
        nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinGeom_;
        }

public:

    //- Runtime type information
    TypeName("mechanicalLaw");


    // Declare run-time constructor selection table

        //- Mechanical law for linear geometry i.e. small strains and rotations
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            linGeomMechLaw,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const nonLinearGeometry::nonLinearType& nonLinGeom
            ),
            (name, mesh, dict, nonLinGeom)
        );

        //- Mechanical law for nonlinear geometry i.e. finite strains and
        //  rotations
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            nonLinGeomMechLaw,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const nonLinearGeometry::nonLinearType& nonLinGeom
            ),
            (name, mesh, dict, nonLinGeom)
        );


    // Selectors

        //- Create a linear geometry mechanical law
        static autoPtr<mechanicalLaw> NewLinGeomMechLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );

        //- Create a nonlinear geometry mechanical law
        static autoPtr<mechanicalLaw> NewNonLinGeomMechLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Constructors

        //- Construct from dictionary
        mechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~mechanicalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return density
        virtual tmp<volScalarField> rho() const = 0;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const = 0;

        //- Return the implicit stiffness surface field
        //  This is the diffusivity for the Laplacian term
        virtual tmp<surfaceScalarField> impKf() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma) = 0;

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update total accumulated fields
        virtual void updateTotalFields()
        {}

        //- Return the desired new time-step
        virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
