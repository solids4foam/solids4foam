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

\*---------------------------------------------------------------------------*/

#include "mechanics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanics, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

mechanics::mechanics
(
    const volTensorField& F,
    const dictionary& dict
)
:
    mesh_(F.mesh()),

    op(mesh_),

    N_(mesh_.Sf()/mesh_.magSf()),

    n_
    (
        IOobject
        (
            "n",
            F.time().timeName(),
            F.db()
        ),
        mesh_.Sf()/mesh_.magSf()
    ),

    S_lm_
    (
        IOobject
        (
            "S_lm",
            F.time().timeName(),
            F.db()
        ),
        F.mesh(),
        dimensionedTensor("S_lm", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    ),

    S_t_
    (
        IOobject
        (
            "S_t",
            F.time().timeName(),
            F.db()
        ),
        F.mesh(),
        dimensionedTensor("S_t", dimensionSet(0,-1,1,0,0,0,0), tensor::zero)
    ),


    cfl_(readScalar(dict.lookup("cfl"))),

    tStep_(0),

    stretch_
    (
        IOobject
        (
            "stretch",
            F.time().timeName(),
            F.db()
        ),
        F.mesh(),
        dimensionedScalar("stretch", dimless, 1.0)
    )

{



    if (cfl_ <= 0.0 || cfl_ > 1.0)
    {
        FatalErrorIn("readControls.H")
            << "Valid type entries are '<= 1' or '> 0' for cfl"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

mechanics::~mechanics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

surfaceVectorField mechanics::spatialNormal(const volTensorField& F_)
{
  surfaceTensorField FcInv = (inv(fvc::interpolate(F_))).ref();
    n_ = (FcInv.T() & N_)/(mag(FcInv.T() & N_));

    return n_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mechanics::correct
(
    const GeometricField<scalar, fvPatchField, volMesh>& Up,
    const GeometricField<scalar, fvPatchField, volMesh>& Us,
    const GeometricField<tensor, fvPatchField, volMesh>& F
)
{
    // Spatial normals
  surfaceTensorField FcInv = (inv(fvc::interpolate(F))).ref();
  surfaceVectorField n_ = ((FcInv.T() & N_)/(mag(FcInv.T() & N_))).ref();

    // Stretch
  volTensorField C_ = (F.T() & F).ref();
    forAll(mesh_.cells(), cell)
    {
        op.eigenStructure(C_[cell]);
        vector eigVal_ = op.eigenValue();
        stretch_[cell] = min( eigVal_.x(), eigVal_.y());
        stretch_[cell] = ::sqrt(min(stretch_[cell], eigVal_.z()));
    }

    if (Pstream::parRun())
    {
        stretch_.correctBoundaryConditions();
    }

    // Stabilisation matrices
    S_lm_ = fvc::interpolate(Up)*n_*n_ + fvc::interpolate(Us)*(I-(n_*n_));

    S_t_ = (n_*n_)/fvc::interpolate(Up) + (I-(n_*n_))/fvc::interpolate(Us);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mechanics::time
(
    Time& runTime,
    dimensionedScalar& deltaT,
    dimensionedScalar Up_time
)
{
    const dimensionedScalar& h = op.minimumEdgeLength();

    deltaT = cfl_*min((h)/Up_time, 0.666*runTime.deltaT());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mechanics::printCentroid() const
{
    vector sum = vector::zero;
    scalar vol = gSum(mesh_.V());

    forAll(mesh_.cells(), cell)
    {
        sum += mesh_.C()[cell]*mesh_.V()[cell];
    }

    if (Pstream::parRun())
    {
        reduce(sum, sumOp<vector>());
    }

    Info << "\nCentroid of geometry = " << sum/vol << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
