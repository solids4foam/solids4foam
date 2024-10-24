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
    coconutModel

Description
    This class ports the models from the CoCoNuT toolbox available
    at https://github.com/pyfsi/coconut.

Author
    Ported to OpenFOAM by Philip Cardiff, UCD.
    See the CoCoNuT toolbox to see the authors or the original Python procedures.

SourceFiles
    coconutModel.C

\*---------------------------------------------------------------------------*/

#ifndef coconutModel_H
#define coconutModel_H

#include "tmp.H"
#include "vectorField.H"
#include "PtrList.H"
#include "Switch.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                        Class coconutModel Declaration
\*---------------------------------------------------------------------------*/

class coconutModel
{
    // Private data

        //- Number of modes to re-use
        const int q_;

        //- Switch to indicate if fields have been added
        Switch added_;

        //- V matrix stored as a list of fields
        PtrList<vectorField> vCurr_;
        PtrList<vectorField> vPrev_;

        //- W matrix stored as a list of fields
        PtrList<vectorField> wCurr_;
        PtrList<vectorField> wPrev_;

        //- r reference values
        vectorField rRef_;

        //- xt reference values
        vectorField xtRef_;


    // Private Member Functions

        //- Unroll a vectorField to a scalarField
        tmp<scalarField> unroll
        (
            const vectorField& x
        ) const;

        //- Roll a scalarField into vectorField
        tmp<vectorField> roll
        (
            const scalarField& x
        ) const;

        //- Convert a PtrList<vectorField> to a scalarRectangularMatrix 
        void convertToMat
        (
            scalarRectangularMatrix& vMat,
            const PtrList<vectorField>& vCurr
        ) const;

        //- Filter modes
        void filter();

        //- Limit
        scalarRectangularMatrix limit
        (
            const scalarRectangularMatrix& matrix,
            const int modes
        ) const;

        //- Combine a and b lists into a matrix and apply the limit function
        void combineAndLimit
        (
            scalarRectangularMatrix& mat,
            const PtrList<vectorField>& aList,
            const PtrList<vectorField>& bList,
            const int modes
        ) const;

        //- Disallow default bitwise copy construct
        coconutModel(const coconutModel&);

        //- Disallow default bitwise assignment
        void operator=(const coconutModel&);


public:

    //- Runtime type information
    TypeName("coconutModel");


    // Constructors

        //- Construct from components
        coconutModel
        (
            const int q    // Number of modes to re-use
        );


    // Destructor

        virtual ~coconutModel()
        {}


    // Member Functions

        // Edit

            //- Initialise values at the start of a time step
            void initializeSolutionStep();

            //- Finalise values at the start of a time step
            void finaliseSolutionStep();

            //- Add
            void add(const vectorField& r, const vectorField& xt);

            //- Is the model ready?
            bool isReady();

            //- Make a prediction
            tmp<vectorField> predict
            (
                const vectorField& r, const int modes = -1
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
