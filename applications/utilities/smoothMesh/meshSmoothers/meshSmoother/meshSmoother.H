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
    newMeshSmoother.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef meshSmoother_H
#define meshSmoother_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

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

        // Mesh smoothing dict
        IOdictionary meshSmootherDict_;

        // Old points instance
        word oldInstance_;

    // Private Member Functions

        //- Disallow copy construct
        meshSmoother(const meshSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const meshSmoother&);


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
                const word name,
                fvMesh& mesh
            ),
            (name, mesh)
        );


    // Selectors

        //- Return a reference to the selected smoother algorithm
        static autoPtr<meshSmoother> New
        (
            const word& name,
                fvMesh& mesh
        );


    // Constructors

        //- Construct from dictionary
        meshSmoother
        (
            const word& name,
                fvMesh& mesh
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

        //- Return reference to meshSmootherDict
        const dictionary& dict() const
        {
            return meshSmootherDict_;
        }

        //- Return reference to mesh
        word& oldInstance()
        {
            return oldInstance_;
        }

        //- The actual smoothing function
        virtual scalar smooth(bool writeIters) = 0;

        //- Mesh writting based on option specified in dictionary
        virtual void writeSmoothedMesh(Time& runTime);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
