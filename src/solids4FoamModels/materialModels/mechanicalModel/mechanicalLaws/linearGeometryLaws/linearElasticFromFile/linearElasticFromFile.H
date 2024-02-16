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

Class
    linearElasticFromFile

Description
    Hookean linear elastic mechanical law.

SourceFiles
    linearElasticFromFile.C

Author
    Philip Cardiff UCD
\*---------------------------------------------------------------------------*/

#ifndef linearElasticFromFile_H
#define linearElasticFromFile_H

#include "mechanicalLaw.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "surfaceFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElasticFromFile Declaration
\*---------------------------------------------------------------------------*/

class linearElasticFromFile
:
    public mechanicalLaw
{
    // Private data

        //- Young's modulus field
        volScalarField E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Shear modulus field
        volScalarField mu_;

        //- 1st Lame parameter field
        volScalarField lambda_;

        //- Shear modulus surface field
        surfaceScalarField muf_;

        //- 1st Lame parameter surface field
        surfaceScalarField lambdaf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticFromFile(const linearElasticFromFile&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticFromFile&);

public:

    //- Runtime type information
    TypeName("linearElasticFromFile");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElasticFromFile
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const label lawI = -1,
            const solidSubMeshes* solidSubMeshes = NULL
        );


    // Destructor

        virtual ~linearElasticFromFile();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> impK() const;

        //- Calculate stress volField
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate stress surfaceField
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
