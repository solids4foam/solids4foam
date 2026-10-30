/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "groundContactContribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(groundContactContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, groundContactContribution, dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


groundContactContribution::groundContactContribution
(
    const word& name,
    const dictionary& dict
)
:
    beamMomentumContribution
    (
        name,
        dict
    ),
    beamMomentumContribDict_(dict.subDict(name + "Coeffs")),
    kNormal_
    (
        readScalar(beamMomentumContribDict_.lookup("kNormal"))
    ),
    cNormal_
    (
        readScalar(beamMomentumContribDict_.lookup("cNormal"))
    ),
    kTangential_
    (
        readScalar(beamMomentumContribDict_.lookup("kTangential"))
    ),
    muFriction_
    (
        readScalar(beamMomentumContribDict_.lookup("muFriction"))
    ),
    groundPlane_(readGroundPlane(beamMomentumContribDict_))
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

plane groundContactContribution::readGroundPlane(const dictionary& dict)
{
    const bool pointAndNormal =
        dict.found("groundPoint") || dict.found("groundNormal");

    if (dict.found("groundZ") == pointAndNormal)
    {
        FatalIOErrorInFunction(dict)
            << "Specify the ground plane either with groundZ, or with"
            << " groundPoint and groundNormal" << exit(FatalIOError);
    }

    if (pointAndNormal)
    {
        // The plane normalises the normal, and stops if it is zero
        return plane
        (
            dict.get<point>("groundPoint"),
            dict.get<vector>("groundNormal")
        );
    }

    // Short form: the horizontal plane z = groundZ
    return plane
    (
        point(0, 0, readScalar(dict.lookup("groundZ"))),
        vector(0, 0, 1)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<scalarSquareMatrix>> groundContactContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<Field<scalarSquareMatrix>> tresult
    (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
    );
    // Field<scalarSquareMatrix>& result = tresult.ref();

    return tresult;
}


tmp<vectorField> groundContactContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Current position of the beam CV cell-centres: the mesh cell-centres
    // plus the reference displacement (registered by the beam model) and the
    // displacement
    const vectorField coord
    (
        mesh.C().internalField()
      + mesh.lookupObject<volVectorField>("refW").internalField()
      + W.internalField()
    );

    // Unit normal of the ground plane, pointing from the ground into the fluid
    const vector& n = groundPlane_.normal();

    // Beam Radius
    const scalar R = bm.R();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    // TO-DO: The friction below only resists sliding along the beam, not
    // sliding across it on the ground

    // Unit tangents to beam centreline at beam CV cell-centres
    const vectorField dRdScellHat(unitTangents(bm));

    // Tangential component of velocity vector
    vectorField Ut
    (
        (
            (U.internalField() & dRdScellHat)
            *dRdScellHat
        )
    );

    vectorField UtHat (Ut/(mag(Ut) + SMALL));

    // Initialise beam cells in contact with ground
    label cellsInContact = 0;

    forAll(result, cellI)
    {
        // Signed distance of the cell-centre from the ground plane,
        // negative below the ground
        const scalar d = groundPlane_.signedDistance(coord[cellI]);

        if (d < 0)
        {
            cellsInContact++;
            // Seabed force on the cell, along the ground normal: the spring
            // pushes the penetrating cell out, and the damper resists its
            // rebound (velocity along the normal). The seabed can only push,
            // so a damper stronger than the spring does not pull the cell in
            const scalar f_gc_normal = max
            (
                (2*R*kNormal_*(-d))
              - (2*R*cNormal_*max(U[cellI] & n, 0)),
                0
            );

            vector f_gc_tangential(vector::zero);

            if
            (
                kTangential_*2.0*R*mag(Ut[cellI])
             >= muFriction_*f_gc_normal
            )
            {
                // Max value of friction as per Coulomb's law
                f_gc_tangential =
                    -muFriction_*f_gc_normal*UtHat[cellI];
            }
            else
            {
                f_gc_tangential = -2.0*R*kTangential_*Ut[cellI];
            }

            // The solver adds external forces to the source with a minus
            // sign, as for gravity, so subtract the seabed force
            result[cellI] -= f_gc_normal*n + f_gc_tangential;
        }
     }

    Info<< "Number of cells in contact : " << cellsInContact << endl;

    return tresult;
}


tmp<vectorField> groundContactContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    // vectorField& result = tresult.ref();

    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
