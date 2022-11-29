/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    differentialStresses

Description
    FunctionObject that outputs the differential stress fields.

    Calculate and write differential stress fields. It assumed that the stress
    tensor is called "sigma" or "sigmaCauchy".

    One scalar fields is created: sigmaDiff
    sigmaDiff is the difference between the most positive/tensile principal stress and the
    most negative/compressive principal stress.

    @verbatim
    differentialStresses1
    {
        // Type of functionObject
        type differentialStresses;

        // Where to load it from (if not already in solver)
        // Note: this is not required if the solver already loads this library.
        functionObjectLibs ("libsolidModelsBekaert.so");

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;
    );
    @endverbatim

Author
    Scott Levie, Philip Cardiff UCD

SourceFiles
    differentialStresses.C

\*---------------------------------------------------------------------------*/

#ifndef differentialStresses_H
#define differentialStresses_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                            Class differentialStresses Declaration
    \*---------------------------------------------------------------------------*/

    class differentialStresses
        : public functionObject
    {
        // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time &runTime_;

        //- Reference to mesh region
        const fvMesh &mesh_;

        //- Reference to compressionPositive flag
        const bool compressionPositive_;

        // Private Member Functions

        //- Calculate and set differential stress vectors
        void calculateEigenValues(
            const symmTensor &sigma);

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        differentialStresses(
            const differentialStresses &);

        //- Disallow default bitwise assignment
        void operator=(const differentialStresses &);

    public:
        //- Runtime type information
        TypeName("differentialStresses");

        // Constructors

        //- Construct from components
        differentialStresses(
            const word &name,
            const Time &,
            const dictionary &);

        // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary &dict);

#ifdef OPENFOAMESIORFOUNDATION
        //- Write
        virtual bool write();
#endif

#ifndef OPENFOAMESIORFOUNDATION
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh &)
        {
        }

        //- Update for changes of mesh
        virtual void movePoints(const pointField &)
        {
        }
#endif
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
