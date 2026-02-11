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
    sphericalCavityAnalyticalSolution

Description
    FunctionObject that generates the analytical displacement and stress fields
    for the spherical cavity (radius = a) in an infinite linear elastic solid
    (Young's modulus of E and Poisson's ration of nu) subjected to a distant
    uniaxial stress T0 applied in the positive Z direction.

    The analytical solution requires the following inputs:
        - farFieldTractionZ : far-field traction in the positive Z direction
        - cavityRadius : radius of the cavity centred on the origin
        - E : Young's modulus
        - nu : Poisson's ratio

    The analytical solution for stress is taken from:

    R.V. Southwell M.A. F.R.S. & H.J. Gough M.B.E. B.Sc. (1926) VI. On the
    concentration of stress in the neighbourhood of a small spherical flaw;
    and on the propagation of fatigue fractures in “Statistically Isotropic”
    materials , The London, Edinburgh, and Dublin Philosophical Magazine and
    Journal of Science: Series 7, 1:1, 71-97, DOI: 10.1080/14786442608633614

    The analytical solution for displacement is taken from:

    Allan F. Bower (2009) Applied Mechanics of Solids, CRC Press, Taylor and
    Francis Group.

    Note: Bower also gives the stress solutions but they seem to have errors
    as the axial (zz) stress formula does not provide the correct stress
    concentration factor at the edge of the cavity.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    sphericalCavityAnalyticalSolution.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalCavityAnalyticalSolution_H
#define sphericalCavityAnalyticalSolution_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sphericalCavityAnalyticalSolution Declaration
\*---------------------------------------------------------------------------*/

class sphericalCavityAnalyticalSolution
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Calculate cell displacement fields
        const Switch cellDisplacement_;

        //- Calculate point displacement fields
        const Switch pointDisplacement_;

        //- Calculate cell stress fields
        const Switch cellStress_;

        //- Calculate point stress fields
        const Switch pointStress_;

        //- Optional history file for writing the errors
        autoPtr<OFstream> historyFilePtr_;

        //- Far field traction
        //  It is assumed that it acts in the positive Z direction
        const scalar T0_;

        //- Cavity radius
        const scalar cavityR_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        sphericalCavityAnalyticalSolution
        (
            const sphericalCavityAnalyticalSolution&
        );

        //- Disallow default bitwise assignment
        void operator=(const sphericalCavityAnalyticalSolution&);


public:

    //- Runtime type information
    TypeName("sphericalCavityAnalyticalSolution");


    // Constructors

        //- Construct from components
        sphericalCavityAnalyticalSolution
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

#ifdef OPENFOAM_NOT_EXTEND
        //- Write
        virtual bool write();
#else
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
