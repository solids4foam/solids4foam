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
    abaqusUmatLinearElastic

Description
    Wrapper class for abaqusUmatLinearElastic.f fortran sub-routine from Abaqus.

SourceFiles
    abaqusUmatLinearElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef abaqusUmatLinearElastic_H
#define abaqusUmatLinearElastic_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class abaqusUmatLinearElastic Declaration
\*---------------------------------------------------------------------------*/

class abaqusUmatLinearElastic
:
    public mechanicalLaw
{
    // Private data

        //- Density
        const dimensionedScalar rho_;

        //- List of material properties
        const List<scalar> properties_;

        //- List of internal state scalar variable fields
        // PtrList<volScalarField> stateVariables_;

        //- "Implicit stiffness" used by the segregated solid models
        //  This only affects the convergence assuming convergence is achieved
        //  For linear elastic solids the ideal value is "2*mu + lambda"
        const dimensionedScalar impK_;

        //- Total strain
        volSymmTensorField epsilon_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        abaqusUmatLinearElastic(const abaqusUmatLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const abaqusUmatLinearElastic&);

public:

    //- Runtime type information
    TypeName("abaqusUmatLinearElastic");


    // Constructors

        //- Construct from dictionary
        abaqusUmatLinearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~abaqusUmatLinearElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return density scale
        virtual const dimensionedScalar& rhoScalar() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar& mu() const;

        //- Bulk modulus
        virtual const dimensionedScalar& K() const;

        //- Young's modulus
        virtual const dimensionedScalar& E() const;

        //- Poisson's ratio
        virtual const dimensionedScalar& nu() const;

        //- Lambda
        virtual const dimensionedScalar& lambda() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
