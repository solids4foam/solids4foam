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
    adaptiveRhieChowStab

Description
    Pressure stabilisation based on the diff-stencil correction operator with
    an adaptive coefficient built from the approximate momentum diagonal.

SourceFiles
    adaptiveRhieChowStab.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveRhieChowStab_H
#define adaptiveRhieChowStab_H

#include "stabilisationModel.H"

namespace Foam
{

class adaptiveRhieChowStab
:
    public stabilisationModel
{
    // Private data

        //- Legacy pressure smoothing coefficient
        const scalar pressureSmoothingCoeff_;

        //- Cached adaptive gamma field
        mutable autoPtr<surfaceScalarField> adaptiveGammaPtr_;

        //- Time index when the adaptive gamma cache was last built
        mutable label adaptiveGammaTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        adaptiveRhieChowStab(const adaptiveRhieChowStab&);

        //- Disallow default bitwise assignment
        void operator=(const adaptiveRhieChowStab&);


public:

    //- Runtime type information
    TypeName("adaptiveRhieChow");


    // Constructors

        //- Construct from dictionary
        adaptiveRhieChowStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~adaptiveRhieChowStab();


    // Member Functions

        //- Re-calculate the scalar stabilisation
        virtual void updateScalar
        (
            const volScalarField& field,
            const volVectorField* gradPtr = nullptr
        ) const;

        //- Provide an implicit approximation of the stabilisation
        virtual const fvScalarMatrix& scalarJacobian
        (
            const volScalarField& field,
            const surfaceScalarField* gammaPtr = nullptr,
            const bool rebuild = false
        ) const;

        //- This model uses an adaptive scalar diffusivity
        virtual bool usesAdaptiveScalarGamma() const
        {
            return true;
        }

        //- Build the adaptive scalar diffusivity from the momentum diagonal
        virtual tmp<surfaceScalarField> adaptiveScalarGamma
        (
            const volVectorField& D,
            const surfaceScalarField& impKf,
            const volScalarField& rho,
            const dimensionedScalar& dampingCoeff
        ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
