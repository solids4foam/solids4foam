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

#include "immersedBody.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "Map.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::triSurface Foam::immersedBody::readSurface
(
    const fvMesh& mesh,
    const fileName& f
)
{
    fileName file(f);
    file.expand();

    // The surface is read by every processor from the case directory
    if (!file.isAbsolute())
    {
        file =
            mesh.time().globalPath()/mesh.time().constant()/"triSurface"/file;
    }

    if (!isFile(file))
    {
        FatalErrorInFunction
            << "Cannot find the immersed body surface " << file
            << exit(FatalError);
    }

    return triSurface(file);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBody::immersedBody
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    mesh_(mesh),
    surface_(readSurface(mesh, dict.get<fileName>("surface"))),
    points0_(surface_.points()),
    searchPtr_(),
    motionPtr_(),
    CofR0_(dict.getOrDefault<point>("CofR", average(points0_))),
    time_(-GREAT),
    internalCells_(),
    surfaceCells_()
{
    // The body is static unless a motion is given
    dictionary motionDict;

    if (dict.found("motion"))
    {
        motionDict = dict.subDict("motion");
    }
    else
    {
        motionDict.add("type", word("static"));
    }

    motionPtr_ = immersedBodyMotion::New(motionDict);

    // The inside test requires a closed surface with outward normals
    if (surface_.nInternalEdges() != surface_.nEdges())
    {
        FatalErrorInFunction
            << "The surface of immersed body " << name_ << " is not closed: "
            << surface_.nEdges() - surface_.nInternalEdges()
            << " edges are not shared by two faces" << exit(FatalError);
    }

    scalar volume = 0;
    forAll(surface_, facei)
    {
        const triPointRef tri(surface_[facei].tri(surface_.points()));
        volume += (tri.centre() & tri.areaNormal())/3;
    }

    if (volume < 0)
    {
        FatalErrorInFunction
            << "The normals of the surface of immersed body " << name_
            << " point into the body: they must point out of it"
            << exit(FatalError);
    }

    Info<< "    Immersed body " << name_ << ": " << surface_.size()
        << " faces, bounding box " << boundBox(points0_, false) << endl;

    move(mesh.time().value());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::immersedBody::CofR() const
{
    return motionPtr_->points(pointField(1, CofR0_), time_)()[0];
}


void Foam::immersedBody::move(const scalar t)
{
    if (searchPtr_ && (!moving() || t == time_))
    {
        return;
    }

    surface_.movePoints(motionPtr_->points(points0_, t));
    searchPtr_.reset(new triSurfaceSearch(surface_));
    time_ = t;
}


void Foam::immersedBody::addOccupancy
(
    volScalarField& lambda,
    const scalar threshold
)
{
    // Cells whose bounding box overlaps the bounding box of the surface
    treeBoundBox surfaceBb(surface_.points());
    surfaceBb.grow(SMALL*mag(surfaceBb.span()));

    const labelList cells(mesh_.cellTree().findBox(surfaceBb));

    // Test each point of these cells once
    const labelListList& cellPoints = mesh_.cellPoints();
    Map<label> pointToIndex(4*cells.size());
    DynamicList<label> points(4*cells.size());

    forAll(cells, i)
    {
        for (const label pointi : cellPoints[cells[i]])
        {
            if (pointToIndex.insert(pointi, points.size()))
            {
                points.append(pointi);
            }
        }
    }

    const boolList pointInside
    (
        searchPtr_->calcInside(pointField(mesh_.points(), points))
    );

    const boolList centreInside
    (
        searchPtr_->calcInside(pointField(mesh_.C(), cells))
    );

    // Occupancy: half the fraction of the cell vertices inside the surface,
    // plus a half if the cell centre is inside the surface
    DynamicList<label> internalCells(cells.size());
    DynamicList<label> surfaceCells(cells.size());
    scalarField& lambdaI = lambda.primitiveFieldRef();

    forAll(cells, i)
    {
        const label celli = cells[i];
        const labelList& curPoints = cellPoints[celli];
        const scalar pointWeight = 0.5/curPoints.size();

        scalar cellLambda = 0;

        for (const label pointi : curPoints)
        {
            if (pointInside[pointToIndex[pointi]])
            {
                cellLambda += pointWeight;
            }
        }

        if (centreInside[i])
        {
            cellLambda += 0.5;
        }

        if (cellLambda > threshold)
        {
            if (cellLambda > 1 - threshold)
            {
                internalCells.append(celli);
            }
            else
            {
                surfaceCells.append(celli);
            }
        }

        lambdaI[celli] = min(max(lambdaI[celli] + cellLambda, 0), 1);
    }

    internalCells_.transfer(internalCells);
    surfaceCells_.transfer(surfaceCells);
}


void Foam::immersedBody::setVelocity(volVectorField& Ui) const
{
    vectorField& UiI = Ui.primitiveFieldRef();

    for (const labelList* cellsPtr : {&internalCells_, &surfaceCells_})
    {
        const labelList& cells = *cellsPtr;

        const vectorField velocity
        (
            motionPtr_->velocity(pointField(mesh_.C(), cells), time_)
        );

        forAll(cells, i)
        {
            UiI[cells[i]] = velocity[i];
        }
    }
}


Foam::vector Foam::immersedBody::force
(
    const volVectorField& f,
    const scalar rho
) const
{
    const scalarField& V = mesh_.V();

    vector F(Zero);

    for (const labelList* cellsPtr : {&internalCells_, &surfaceCells_})
    {
        for (const label celli : *cellsPtr)
        {
            F -= f[celli]*V[celli];
        }
    }

    reduce(F, sumOp<vector>());

    return rho*F;
}


Foam::vector Foam::immersedBody::torque
(
    const volVectorField& f,
    const scalar rho
) const
{
    const scalarField& V = mesh_.V();
    const vectorField& C = mesh_.C();
    const point CofR(this->CofR());

    vector T(Zero);

    for (const labelList* cellsPtr : {&internalCells_, &surfaceCells_})
    {
        for (const label celli : *cellsPtr)
        {
            T -= ((C[celli] - CofR) ^ f[celli])*V[celli];
        }
    }

    reduce(T, sumOp<vector>());

    return rho*T;
}


Foam::vector Foam::immersedBody::momentum
(
    const volVectorField& U,
    const volScalarField& lambda
) const
{
    const scalarField& V = mesh_.V();

    vector M(Zero);

    for (const labelList* cellsPtr : {&internalCells_, &surfaceCells_})
    {
        for (const label celli : *cellsPtr)
        {
            M += lambda[celli]*U[celli]*V[celli];
        }
    }

    reduce(M, sumOp<vector>());

    return M;
}


void Foam::immersedBody::writeSurface(const fileName& file) const
{
    if (Pstream::master())
    {
        mkDir(file.path());
        surface_.write(file);
    }
}


// ************************************************************************* //
