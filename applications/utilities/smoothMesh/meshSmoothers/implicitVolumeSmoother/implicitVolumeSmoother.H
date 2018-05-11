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
    implicitVolumeSmoother

Description
    Volume smoothing, where each point is moved based on the neighbouring cell
    centres using weights based on the volumes of the cells.

    The smoothing is performed implicitly by solving a system of nonlinear
    equations.

    The volFields can be advected using an explicit advection transport equation
    where the mesh Courant number can be limited.

    Current issue for parallel in 2-D: a processor may not have cells in the
    subSet, causing it to give an error when calculating the empty directions.

SourceFiles
    implicitVolumeSmoother.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef implicitVolumeSmoother_H
#define implicitVolumeSmoother_H

#include "meshSmoother.H"
#include "newFvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class implicitVolumeSmoother Declaration
\*---------------------------------------------------------------------------*/

class implicitVolumeSmoother
:
    public meshSmoother
{
    // Private Data

        //- Sub-mesh containing the cells to be smoothed
        mutable autoPtr<newFvMeshSubset> subMeshToSmooth_;


    // Private Member Functions

        //- Return a reference to the sub mesh to be smoothed
        fvMesh& subMeshToSmooth();

        //- Return a const reference to the sub mesh to be smoothed
        const fvMesh& subMeshToSmooth() const;

        //- Return a reference to the subset mesh to be smoothed
        newFvMeshSubset& subsetMeshToSmooth();

        //- Return a const reference to the subset mesh to be smoothed
        const newFvMeshSubset& subsetMeshToSmooth() const;

        //- Return a reference to the sub mesh to be smoothed
        void makeSubMeshToSmooth() const;

        //- Calculate and return the mesh diffusivity field
        tmp<surfaceScalarField> diffusivity() const;

        //- Calculate and return the cell centres field with appropriate
        //  boundary conditions
        tmp<volVectorField> cellCentresField() const;

        //- Calculate and return the cell volumes field with appropriate
        //  boundary conditions
        tmp<volScalarField> cellVolumesField() const;

        //- Map motionD from the cell centres to the vertices
        void mapVolToPoint
        (
            const volVectorField& cellMotionD,
            pointVectorField& pointMotionD
        );

        //- Calculate mesh point motion using implicit cell-based method
        void calculatePointMotionImplicitly(pointVectorField& pointMotionD);

        //- Disallow default bitwise copy construct
        implicitVolumeSmoother(const implicitVolumeSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const implicitVolumeSmoother&);

public:

    //- Runtime type information
    TypeName("implicitVolume");

    // Static data members


    // Constructors

        //- Construct from dictionary
        implicitVolumeSmoother
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~implicitVolumeSmoother();


    // Member Functions

        //- Clear out demand driven data
        virtual void clearOut();

        //- Smoothing function
        virtual scalar smooth();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
