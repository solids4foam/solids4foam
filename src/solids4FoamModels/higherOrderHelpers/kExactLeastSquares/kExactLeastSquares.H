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
    Foam::kExactLeastSquares


Description
    The class constructs local polynomial approximations of a scalar field
    whose stored cell values represent cell volume averages. A weighted
    k-exact least-squares procedure is applied over each cell stencil to
    reconstruct a polynomial of order N about the owner-cell centre.
    Neighbour observations are formed from differences between cell averages.
    Consequently, a polynomial of degree up to N is reproduced
    exactly when the reconstruction system has sufficient rank.

    The resulting coefficients allow reconstruction of the field value,
    gradient, and second and third derivatives. Reconstruction coefficients
    are precomputed and stored for reuse during the simulation, allowing
    efficient evaluation at cell centres, arbitrary points, and face
    quadrature points during the discretisation phase.

    For each cell, the class solves a small dense weighted least-squares system
    in which the matrix contains volume-averaged, mean-free polynomial basis
    functions evaluated over the stencil cells. At internal faces, gradients
    from the owner and neighbour reconstructions are averaged at each
    quadrature point. Prescribed-value boundaries use the boundary values at
    the quadrature points as additional weighted observations, while symmetry
    boundaries use reflected cell observations.

    Sources:
        Sejekan CB, Ollivier-Gooch CF. Improving finite-volume diffusive
        fluxes through better reconstruction. Computers & Fluids, 139:
        216-232, 2016.

        Ollivier-Gooch C, Van Altena M. A High-Order-Accurate Unstructured
        Mesh Finite-Volume Scheme for the Advection-Diffusion Equation.
        Journal of Computational Physics, 181:729-752, 2002.

SourceFiles
    kExactLeastSquares.C
    kExactLeastSquaresTemplates.C

    The per-cell polynomial machinery is implemented in cellLeastSquares and
    shared with leastSquares, the point-value scheme; this class selects the
    cell-average representation through cellAverageUnknowns().

SourceFiles
    kExactLeastSquares.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef kExactLeastSquares_H
#define kExactLeastSquares_H

#include "cellLeastSquares.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class kExactLeastSquares Declaration
\*---------------------------------------------------------------------------*/

class kExactLeastSquares
:
    public cellLeastSquares
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        kExactLeastSquares(const kExactLeastSquares&);

        //- Disallow default bitwise assignment
        void operator=(const kExactLeastSquares&);


public:

    //- Runtime type information
    TypeName("kExactLeastSquares");


    // Constructors

        //- Construct from components
        kExactLeastSquares
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );


    //- Destructor
    virtual ~kExactLeastSquares();


    // Member Functions

        //- The stored cell values are cell averages
        virtual bool cellAverageUnknowns() const
        {
            return true;
        }

        //- Return second-order central cell moments owned by quadrature
        const List<symmTensor>& secondOrderCellMoments() const;

        //- Return third-order central cell moments owned by quadrature
        const List<symmTensor3rdOrder>& thirdOrderCellMoments() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
