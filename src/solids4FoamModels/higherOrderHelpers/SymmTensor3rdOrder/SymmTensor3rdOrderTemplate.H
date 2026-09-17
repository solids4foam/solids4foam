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
    Foam::SymmTensor3rdOrder

Description
    Templated fully symmetric third-order tensor with 10 independent
    components in three-dimensional space.

    The class represents a rank-3 tensor that is symmetric with respect
    to permutation of all indices, i.e.

        T_ijk = T_jik = T_ikj = T_kij = T_jki = T_kji

    and therefore stores only the following independent components:

        xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz

    The class is intended for use in high-order finite-volume and moving
    least squares (MLS) reconstruction procedures, where it naturally
    appears as the tensor of third derivatives of a scalar field or as
    the coefficient tensor of the cubic term in a local polynomial
    expansion.

SourceFiles
    SymmTensor3rdOrderTemplateI.H

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef SymmTensor3rdOrderTemplate_H
#define SymmTensor3rdOrderTemplate_H

#include "VectorSpace.H"
#ifdef FOAMEXTEND
    #include "symmTensor.H"
#else
    #include "SymmTensor.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SymmTensor3rdOrder Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SymmTensor3rdOrder
:
    public VectorSpace<SymmTensor3rdOrder<Cmpt>, Cmpt, 10>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef SymmTensor3rdOrder<label> labelType;


    // Member constants

        enum
        {
            rank = 3
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

        static const SymmTensor3rdOrder zero;
        static const SymmTensor3rdOrder one;
        static const SymmTensor3rdOrder max;
        static const SymmTensor3rdOrder min;


        //- Component labeling enumeration
        enum components
        {
            XXX,
            XXY,
            XXZ,
            XYY,
            XYZ,
            XZZ,
            YYY,
            YYZ,
            YZZ,
            ZZZ
        };


    // Constructors

        //- Construct null
        inline SymmTensor3rdOrder();

        //- Construct given VectorSpace
        inline SymmTensor3rdOrder
        (
            const VectorSpace<SymmTensor3rdOrder<Cmpt>, Cmpt, 10>&
        );

        //- Construct given SymmTensor (for compilation)
        //inline SymmTensor3rdOrder(const SymmTensor<Cmpt>&);

        //- Construct given the ten components
        inline SymmTensor3rdOrder
        (
            const Cmpt txxx,
            const Cmpt txxy,
            const Cmpt txxz,
            const Cmpt txyy,
            const Cmpt txyz,
            const Cmpt txzz,
            const Cmpt tyyy,
            const Cmpt tyyz,
            const Cmpt tyzz,
            const Cmpt tzzz
        );

        //- Construct from Istream
        inline SymmTensor3rdOrder(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xxx() const;
            inline const Cmpt& xxy() const;
            inline const Cmpt& xxz() const;
            inline const Cmpt& xyy() const;
            inline const Cmpt& xyz() const;
            inline const Cmpt& xzz() const;
            inline const Cmpt& yyy() const;
            inline const Cmpt& yyz() const;
            inline const Cmpt& yzz() const;
            inline const Cmpt& zzz() const;

            inline Cmpt& xxx();
            inline Cmpt& xxy();
            inline Cmpt& xxz();
            inline Cmpt& xyy();
            inline Cmpt& xyz();
            inline Cmpt& xzz();
            inline Cmpt& yyy();
            inline Cmpt& yyz();
            inline Cmpt& yzz();
            inline Cmpt& zzz();


        //- Assign to a SymmTensor
        //- maybe I need this to compile
        //inline void operator=(const SymmTensor<Cmpt>&);


};

template<class Cmpt>
inline Cmpt cubicForm
(
    const SymmTensor3rdOrder<Cmpt>& T,
    const Vector<Cmpt>& d
);

template<class Cmpt>
class symmTypeOfRank<Cmpt, 3>
{
public:

    typedef SymmTensor3rdOrder<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SymmTensor3rdOrderTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
