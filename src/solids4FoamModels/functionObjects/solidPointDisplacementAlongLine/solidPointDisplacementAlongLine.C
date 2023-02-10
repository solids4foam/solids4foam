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

\*----------------------------------------------------------------------------*/

#include "solidPointDisplacementAlongLine.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "volPointInterpolation.H"
#else
    #include "newLeastSquaresVolPointInterpolation.H"
#endif
#ifdef OPENFOAMFOUNDATION
    #include "OSspecific.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPointDisplacementAlongLine, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPointDisplacementAlongLine,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPointDisplacementAlongLine::writeData()
{
    
    // Lookup the solid mesh
    const fvMesh& mesh = time_.lookupObject<fvMesh>(region_);

    if (mesh.foundObject<pointVectorField>("pointD"))
    {
        // Lookup the point displacement field
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        vector pointDValue = vector::zero;
        
        //Obtain pointD for all values on the specified line
        for (scalar i = 0; i < no_points_; i++)
        {
			if (pointIDSorted_[i] > -1)
			{
				pointDValue = pointD[pointIDSorted_[i]];
			}
			reduce(pointDValue, sumOp<vector>());

			if (Pstream::master())
			{
				historyFilePtr_()
				<< pointIDSorted_[i] << " " << pointCoordSorted_[i]
				<< " " << pointDValue.x()
				<< " " << pointDValue.y()
				<< " " << pointDValue.z()
				<< " " << mag(pointDValue)
				<< endl;
			}
	    }
    }
    else
    {
        InfoIn(this->name() + " function object constructor")
            << "pointVectorField pointD not found" << endl;
    }

    return true;
}

auto Foam::solidPointDisplacementAlongLine::sortByCompX(DynamicList<vector> pointCoord_, DynamicList<label> pointID_)
{
	for (scalar i = 0; i < pointCoord_.size(); i++)
	{
		for (scalar k = i + 1; k < pointCoord_.size(); k++)
		{
			if (mag(pointCoord_[i].component(vector::X)) > mag(pointCoord_[k].component(vector::X)))
			{
				vector tempCoord(pointCoord_[i]);
				pointCoord_[i] = pointCoord_[k];
				pointCoord_[k] = tempCoord;
				
				label tempID(pointID_[i]);
				pointID_[i] = pointID_[k];
				pointID_[k] = tempID;
			} 
		}
	}
	
	return Tuple2<DynamicList<vector>, DynamicList<label>>(pointCoord_, pointID_);
}

auto Foam::solidPointDisplacementAlongLine::sortByCompY(DynamicList<vector> pointCoord_, DynamicList<label> pointID_)
{
	for (scalar i = 0; i < pointCoord_.size(); i++)
	{
		for (scalar k = i + 1; k < pointCoord_.size(); k++)
		{
			if (mag(pointCoord_[i].component(vector::Y)) > mag(pointCoord_[k].component(vector::Y)))
			{
				vector tempCoord(pointCoord_[i]);
				pointCoord_[i] = pointCoord_[k];
				pointCoord_[k] = tempCoord;
				
				label tempID(pointID_[i]);
				pointID_[i] = pointID_[k];
				pointID_[k] = tempID;
			} 
		}
	}
	
	return Tuple2<DynamicList<vector>, DynamicList<label>>(pointCoord_, pointID_);
}

auto Foam::solidPointDisplacementAlongLine::sortByCompZ(DynamicList<vector> pointCoord_, DynamicList<label> pointID_)
{
	for (scalar i = 0; i < pointCoord_.size(); i++)
	{
		for (scalar k = i + 1; k < pointCoord_.size(); k++)
		{
			if (mag(pointCoord_[i].component(vector::Z)) > mag(pointCoord_[k].component(vector::Z)))
			{
				vector tempCoord(pointCoord_[i]);
				pointCoord_[i] = pointCoord_[k];
				pointCoord_[k] = tempCoord;
				
				label tempID(pointID_[i]);
				pointID_[i] = pointID_[k];
				pointID_[k] = tempID;
			
			} 
		}
	}
	
	return Tuple2<DynamicList<vector>, DynamicList<label>>(pointCoord_, pointID_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPointDisplacementAlongLine::solidPointDisplacementAlongLine
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    region_(dict.lookupOrDefault<word>("region", "UNDEFINED")),
    pointID_(-1),
    no_points_(0),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    // Set region if it is undefined
    if (region_ == "UNDEFINED")
    {
        if (time_.foundObject<fvMesh>("solid"))
        {
            region_ = "solid";
        }
        else
        {
            region_ = fvMesh::defaultRegion;
        }
    }
    Info<< "    region = " << region_ << endl;

    // Lookup the start and end point of the line
    const vector pointA(dict.lookup("StartPoint"));
    const vector pointB(dict.lookup("EndPoint"));

    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // Find the closest point on the line
        scalar minDist = 1e-6;
        //Define vector between points A and B
        const vector line_vector(pointB.x() - pointA.x(), pointB.y() - pointA.y(), pointB.z() - pointA.z());
        
        forAll(mesh.points(), pI)
        {
			//Define vector between point A and mesh point P
			vector point_vector(mesh.points()[pI].component(vector::X) - pointA.x(), mesh.points()[pI].component(vector::Y) - pointA.y(), mesh.points()
			[pI].component(vector::Z) - pointA.z());
    	    
    	    //Check whether point is within the region defined by the segment
	        if ( mag(mesh.points()[pI].component(vector::X)) <= mag(pointB.x()) && mag(mesh.points()[pI].component(vector::X)) >= mag(pointA.x()) && 
                 mag(mesh.points()[pI].component(vector::Y)) <= mag(pointB.y()) && mag(mesh.points()[pI].component(vector::Y)) >= mag(pointA.y()) &&
                 mag(mesh.points()[pI].component(vector::Z)) <= mag(pointB.z()) && mag(mesh.points()[pI].component(vector::Z)) >= mag(pointA.z()) )
            {   
			   //Calculate coordinates of projection point on line
			   vector proj_pt(((point_vector & line_vector)/mag(line_vector))*(line_vector/mag(line_vector)) + pointA);
               
               //Calculate distance between mesh point and projection point
               scalar dist = mag(mesh.points()[pI] - proj_pt);
               
               //Check if mesh point is on the line
               if (dist < minDist)
               {
                   pointID_.append(pI);
                   pointCoord_.append(mesh.points()[pI]); 
                   no_points_++;
               }
            }

		    // Find global closest point
		    const scalar globalMinDist = returnReduce(minDist, minOp<scalar>());
		    int procNo = -1;
		    if (mag(globalMinDist - minDist) < SMALL)
		    {
		        procNo = Pstream::myProcNo();
		    }
		    else
		    {
		        pI = -1;
		    }

		    // More than one processor can have the point so we will take the proc
		    // with the highest processor number
		    const int globalMaxProc = returnReduce(procNo, maxOp<int>());
		    if (mag(globalMaxProc - procNo) > SMALL)
		    {
		        pI = -1;
		    }

		    if (pI > -1)
		    {
		        Pout<< "    distance from specified point is " << minDist
		            << endl;
		    }

		    if (returnReduce(pI, maxOp<int>()) == -1)
		    {
		        FatalErrorIn("solidPointDisplacementAlongLine::solidPointDisplacementAlongLine")
		            << "Something went wrong: no proc found a point!"
		            << abort(FatalError);
		    }
	    
	    }
	    
	    //Sort point coordinates by x, y or z-coordinates
	    DynamicList<vector> pointCoordSorted();
	    DynamicList<label> pointIDSorted();
	    
	    if (mag(pointCoord_[0].component(vector::Y)) == mag(pointCoord_[1].component(vector::Y)) && mag(pointCoord_[0].component(vector::Z)) ==
	    mag(pointCoord_[1].component(vector::Z))) 
	    {
			//Sort point coordinates by x-coordinates
			pointCoordSorted_=sortByCompX(pointCoord_, pointID_).first();
			pointIDSorted_=sortByCompX(pointCoord_, pointID_).second();
			
			Info << "SortByCompX called" << pointIDSorted_.size() << pointCoordSorted_.size() << endl;
	    }
	    else if (mag(pointCoord_[0].component(vector::X)) == mag(pointCoord_[1].component(vector::X)) && mag(pointCoord_[0].component(vector::Z)) ==
	    mag(pointCoord_[1].component(vector::Z)))
	    {
	    	//Sort point coordinates by y-coordinates
			pointCoordSorted_=sortByCompY(pointCoord_, pointID_).first();
			pointIDSorted_=sortByCompY(pointCoord_, pointID_).second();
			
			Info << "SortByCompY called" << endl;
	    }
	    else if (mag(pointCoord_[0].component(vector::X)) == mag(pointCoord_[1].component(vector::X)) && mag(pointCoord_[0].component(vector::Y)) ==
	    mag(pointCoord_[1].component(vector::Y)))
	    {
	    	//Sort point coordinates by z-coordinates
			pointCoordSorted_=sortByCompZ(pointCoord_, pointID_).first();
			pointIDSorted_=sortByCompZ(pointCoord_, pointID_).second();
			
			Info << "SortByCompZ called" << endl;
	    }
	    else 
	    {
	    	//Sort point coordinates by x-coordinates
			pointCoordSorted_=sortByCompX(pointCoord_, pointID_).first();
			pointIDSorted_=sortByCompX(pointCoord_, pointID_).second();
			Info << "SortByCompX since two or three coordinates are different" << endl;
	    }
	    
	    // File update
	    if (Pstream::master())
	    {
	        fileName historyDir;

	        const word startTimeName =
	            time_.timeName(mesh.time().startTime().value());

	        if (Pstream::parRun())
	        {
	            // Put in undecomposed case (Note: gives problems for
	            // distributed data running)
	            historyDir = time_.path()/".."/"postProcessing"/startTimeName;
	        }
	        else
	        {
	            historyDir = time_.path()/"postProcessing"/startTimeName;
	        }

	        // Create directory if does not exist.
	        mkDir(historyDir);

	        // Open new file at start up
	        historyFilePtr_.reset
	        (
	            new OFstream
	            (
	                historyDir/"solidPointDisplacementAlongLine_" + name + ".dat"
	            )
	        );

	        // Add headers to output data
	        if (historyFilePtr_.valid())
	        {
	            historyFilePtr_()
	                << "# PointID PointCoord Dx Dy Dz magD" << endl;
	        }
	    }
    }    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPointDisplacementAlongLine::start()
{
    return false;
}

#if FOAMEXTEND
bool Foam::solidPointDisplacementAlongLine::execute(const bool forceWrite)
#else
bool Foam::solidPointDisplacementAlongLine::execute()
#endif
{
    return writeData();
}


bool Foam::solidPointDisplacementAlongLine::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::solidPointDisplacementAlongLine::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
