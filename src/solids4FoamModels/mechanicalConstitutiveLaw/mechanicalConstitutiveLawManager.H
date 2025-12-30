/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mechanicalConstitutiveLawManager

Description
    Manager class responsible for applying mechanical constitutive laws
    to the computational fields.

    The mechanicalConstitutiveLawManager owns one or more mechanical
    constitutive laws together with their associated constitutive state
    objects. It is responsible for extracting kinematic fields from the
    solver, constructing the appropriate kinematics containers, and
    invoking constitutive evaluation.

    This class forms the interface between OpenFOAM field types
    (e.g. volTensorField, pointTensorField) and the discretisation-
    agnostic mechanicalConstitutiveLaw interface. All concerns related
    to field location, mesh topology, and region selection are handled
    exclusively at this level.

    In its initial form, the manager supports a single constitutive law
    applied to the entire domain. This functionality may be extended
    incrementally to support multiple laws, multiple regions, boundary
    contributions, and tangent-related quantities.

    Notes:
    - Constitutive laws are stateless; all history variables are stored
      in mechanicalConstitutiveLawState objects owned by the manager.
    - The manager is the only component that depends on OpenFOAM
      GeometricField types.

SourceFiles
    mechanicalConstitutiveLawManager.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawManager_H
#define mechanicalConstitutiveLawManager_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class mechanicalConstitutiveLawManager Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawManager
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Constitutive laws
        PtrList<mechanicalConstitutiveLaw> laws_;

        //- Constitutive states
        PtrList<mechanicalConstitutiveLawState> states_;

        //- Constitutive state per law per boundary patch
        //  boundaryStates_[lawI][patchI] has size equal to the number of
        //  boundary faces on patch patchI belonging to law lawI
        List<PtrList<mechanicalConstitutiveLawState>> boundaryStates_;

        //- Cell indices per mechanical constitutive law
        //  Each entry corresponds to one mechanicalConstitutiveLaw and contains
        //  the list of local cell indices to which that law applies. The union of
        //  all lists must cover all cells in the mesh, and no cell may appear in
        //  more than one list.
        List<labelList> lawCells_;

        //- Boundary patch face indices per mechanical constitutive law.
        //
        //  For each mechanicalConstitutiveLaw, this stores a list over boundary
        //  patches, where each entry contains the local face indices (within that
        //  patch) whose owner cell belongs to the corresponding law.
        //
        //  In other words:
        //    lawBoundaryFaces_[lawI][patchI] gives the list of boundary faces on
        //    patch patchI whose adjacent cell is assigned to law lawI.
        List<List<labelList>> lawBoundaryFaces_;


    // Private Member Functions

        //- No copy construct
        mechanicalConstitutiveLawManager
        (
            const mechanicalConstitutiveLawManager&
        ) = delete;

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawManager&) = delete;


public:


    // Constructors

        //- Construct from a mesh and dict
        mechanicalConstitutiveLawManager
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~mechanicalConstitutiveLawManager();


    // Member Functions

    // Edit

        //- Update the Cauchy stress field using the configured mechanical
        //  constitutive laws for a small-strain formulation.
        //  The stress is computed from the current and previous displacement
        //  gradients and the supplied time increment. The stress field is
        //  updated in-place. Constitutive history variables are stored and
        //  updated internally via mechanicalConstitutiveLawState objects.
        void updateStress
        (
            const volTensorField& gradD,
            const volTensorField& gradD0,
            const scalar dt,
            volSymmTensorField& stress
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
