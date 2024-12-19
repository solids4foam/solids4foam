/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::angularMomentum

Description
    Angular momentum preservation algorithm

SourceFiles
    angularMomentum.C

\*---------------------------------------------------------------------------*/

#ifndef angularMomentum_H
#define angularMomentum_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class angularMomentum
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Material density
        dimensionedScalar rho_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        angularMomentum(const angularMomentum&);

        //- Disallow default bitwise assignment
        void operator=(const angularMomentum&);

public:

    //- Runtime type information
    TypeName("angularMomentum");


    // Constructors

        //- Construct from mesh and dictionary
        explicit angularMomentum
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~angularMomentum();


    // Member functions

        // Edit

            //- Angular momentum conservation
            void AMconservation
            (
                GeometricField<vector, fvPatchField, volMesh>&,
                GeometricField<vector, fvPatchField, volMesh>&,
                const GeometricField<vector, fvPatchField, volMesh>&,
                const scalar&,
                const dimensionedScalar& pDeltaT 
            ) const;

        // Print

            //- Print global momentum
            void printGlobalMomentum
            (
                const GeometricField<vector, fvPatchField, volMesh>&,
                const GeometricField<vector, fvPatchField, volMesh>&
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //