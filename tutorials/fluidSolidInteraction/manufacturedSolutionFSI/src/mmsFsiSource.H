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
    Foam::fv::mmsFsiSource

Description
    Momentum source of the manufactured FSI solution.

    - region fluid: the source per unit mass of the kinematic momentum
      equation, evaluated at the cell centres of the current (moving) fluid
      mesh;
    - region solid: the body force per unit reference volume of the total
      Lagrangian momentum equation, evaluated at the cell centres of the
      undeformed solid mesh.

    The source is explicit. It is re-evaluated when the time or the cell
    centres change.

Usage
    \verbatim
    mmsSource
    {
        type        mmsFsiSource;
        region      fluid;    // or solid
        fields      (U);      // (D) for the solid
    }
    \endverbatim

SourceFiles
    mmsFsiSource.C

\*---------------------------------------------------------------------------*/

#ifndef mmsFsiSource_H
#define mmsFsiSource_H

#include "fvOption.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class mmsFsiSource Declaration
\*---------------------------------------------------------------------------*/

class mmsFsiSource
:
    public fv::option
{
    // Private Data

        //- Region of the source: fluid or solid
        const word region_;

        //- Time at which the source was last evaluated
        mutable scalar sourceTime_;

        //- Cell centres at which the source was last evaluated
        mutable vectorField sourceC_;

        //- Source at the cell centres
        mutable vectorField source_;


    // Private Member Functions

        //- Return the source, evaluated again if the time or the cell
        //  centres changed
        const vectorField& source() const;


public:

    //- Runtime type information
    TypeName("mmsFsiSource");


    // Constructors

        //- Construct from explicit source name and mesh
        mmsFsiSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        mmsFsiSource(const mmsFsiSource&) = delete;

        //- No copy assignment
        void operator=(const mmsFsiSource&) = delete;


    //- Destructor
    virtual ~mmsFsiSource() = default;


    // Member Functions

        //- Add the explicit source to the momentum equation
        virtual void addSup(fvMatrix<vector>& eqn, const label fieldi);

        //- Add the explicit source to the density-weighted momentum
        //  equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
