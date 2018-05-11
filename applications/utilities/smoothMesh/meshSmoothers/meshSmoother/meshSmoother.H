/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    meshSmoother

Description
    Smoothing mesh

SourceFiles
    meshSmoother.C
    meshSmootherTemplates.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef meshSmoother_H
#define meshSmoother_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
//#include "conservativeMeshToMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//    class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class meshSmoother Declaration
\*---------------------------------------------------------------------------*/

class meshSmoother
{
    // Private data

        // Reference to the mesh
        fvMesh& mesh_;

        // Take a copy of the meshSmooth dict
        const dictionary dict_;

        // Old points instance
        word oldInstance_;

        //- Store a copy of the old mesh points
        //  This is used by the mapFields function
        pointField oldPoints_;

        //- Flag to record if the advectFields function has been called
        bool fieldsAdvected_;


    // Private Member Functions

        //- Map vol fields
        // template<class Type>
        // void MapConservativeVolFields
        // (
        //     const conservativeMeshToMesh& meshToMeshInterp,
        //     const label method
        // );

        //- Advect vol fields
        template<class Type>
        void AdvectVolFields
        (
            const surfaceScalarField& sweptVol,
            const volScalarField& volOld,
            const volScalarField& volNew,
            const wordList& ignoreFields
        );

        //- Disallow copy construct
        meshSmoother(const meshSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const meshSmoother&);


protected:

        //- Correct boundary motion of given point motion field
        void correctBoundaryMotion(pointVectorField& pointMotionD);


public:

    //- Runtime type information
    TypeName("meshSmoother");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            meshSmoother,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Return a pointer to the new smoother
        static autoPtr<meshSmoother> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from mesh and dictionary
        meshSmoother
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~meshSmoother()
        {}


    // Member Functions

        //- Return reference to mesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to meshSmootherDict
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return reference to mesh
        word& oldInstance()
        {
            return oldInstance_;
        }

        //- Should the mesh and field be overwritten
        const Switch overwrite() const
        {
            return Switch(dict_.lookupOrDefault<Switch>("overwrite", false));
        }

        //- Clear out demand driven data
        virtual void clearOut()
        {}

        //- Fields advected
        bool fieldsAdvected() const
        {
            return fieldsAdvected_;
        }

        //- The actual smoothing function
        virtual scalar smooth() = 0;

        //- Store the old points before the smoothing is performed
        virtual void storeOldPoints(const pointField& oldPoints);

        //- Map all GeometricFields (e.g volFields, pointFields, surfaceFields)
        //  associated with the mesh
        virtual void mapFields()
        {
            notImplemented(type() + "::mapFields()");
        }

        //- Advect all the GeometricFields (e.g volFields, pointFields,
        //  surfaceFields) using the suplied mesh flux
        virtual void advectFields
        (
            const scalarField& sweptVol,
            const scalarField& volOld,
            const scalarField& volNew
        );

        //- Write the mesh
        virtual void writeMesh();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshSmootherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
