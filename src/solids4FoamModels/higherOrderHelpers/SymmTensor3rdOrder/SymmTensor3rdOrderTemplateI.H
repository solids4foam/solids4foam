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

#include "SymmTensor.H"
#include "Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<class Cmpt>
inline SymmTensor3rdOrder<Cmpt>::SymmTensor3rdOrder()
{
    this->v_[XXX] = 0;
    this->v_[XXY] = 0;
    this->v_[XXZ] = 0;
    this->v_[XYY] = 0;
    this->v_[XYZ] = 0;
    this->v_[XZZ] = 0;
    this->v_[YYY] = 0;
    this->v_[YYZ] = 0;
    this->v_[YZZ] = 0;
    this->v_[ZZZ] = 0;
}


// Construct given VectorSpace
template<class Cmpt>
inline SymmTensor3rdOrder<Cmpt>::SymmTensor3rdOrder
(
    const VectorSpace<SymmTensor3rdOrder<Cmpt>, Cmpt, 10>& vs
)
:
    VectorSpace<SymmTensor3rdOrder<Cmpt>, Cmpt, 10>(vs)
{}


// Construct given the ten Cmpts
template<class Cmpt>
inline SymmTensor3rdOrder<Cmpt>::SymmTensor3rdOrder
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
)
{
    this->v_[XXX] = txxx;
    this->v_[XXY] = txxy;
    this->v_[XXZ] = txxz;
    this->v_[XYY] = txyy;
    this->v_[XYZ] = txyz;
    this->v_[XZZ] = txzz;
    this->v_[YYY] = tyyy;
    this->v_[YYZ] = tyyz;
    this->v_[YZZ] = tyzz;
    this->v_[ZZZ] = tzzz;
}


// Construct from Istream
template<class Cmpt>
inline SymmTensor3rdOrder<Cmpt>::SymmTensor3rdOrder(Istream& is)
:
    VectorSpace<SymmTensor3rdOrder<Cmpt>, Cmpt, 10>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::xxx() const
{
    return this->v_[XXX];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::xxy() const
{
    return this->v_[XXY];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::xxz() const
{
    return this->v_[XXZ];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::xyy() const
{
    return this->v_[XYY];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::xyz() const
{
    return this->v_[XYZ];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::xzz() const
{
    return this->v_[XZZ];
}


template<class Cmpt>
inline const Cmpt&  SymmTensor3rdOrder<Cmpt>::yyy() const
{
    return this->v_[YYY];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::yyz() const
{
    return this->v_[YYZ];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::yzz() const
{
    return this->v_[YZZ];
}


template<class Cmpt>
inline const Cmpt& SymmTensor3rdOrder<Cmpt>::zzz() const
{
    return this->v_[ZZZ];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::xxx()
{
    return this->v_[XXX];
}

template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::xxy()
{
    return this->v_[XXY];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::xxz()
{
    return this->v_[XXZ];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::xyy()
{
    return this->v_[XYY];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::xyz()
{
    return this->v_[XYZ];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::xzz()
{
    return this->v_[XZZ];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::yyy()
{
    return this->v_[YYY];
}


template<class Cmpt>
inline Cmpt&  SymmTensor3rdOrder<Cmpt>::yyz()
{
    return this->v_[YYZ];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::yzz()
{
    return this->v_[YZZ];
}


template<class Cmpt>
inline Cmpt& SymmTensor3rdOrder<Cmpt>::zzz()
{
    return this->v_[ZZZ];
}


template<class Cmpt>
inline Cmpt cubicForm
(
    const SymmTensor3rdOrder<Cmpt>& T,
    const Vector<Cmpt>& d
)
{
    const scalar dx = d.x();
    const scalar dy = d.y();
    const scalar dz = d.z();

    return
        T[0]*pow3(dx)
      + 3*T[1]*dx*dx*dy
      + 3*T[2]*dx*dx*dz
      + 3*T[3]*dx*dy*dy
      + 6*T[4]*dx*dy*dz
      + 3*T[5]*dx*dz*dz
      +   T[6]*pow3(dy)
      + 3*T[7]*dy*dy*dz
      + 3*T[8]*dy*dz*dz
      +   T[9]*pow3(dz);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Dot product between 3rd order tensor and vector returns 2nd order symmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&(const SymmTensor3rdOrder<Cmpt>& t, const Vector<Cmpt>& v)
{
    return SymmTensor<Cmpt>
    (
        t.xxx()*v.x() + t.xxy()*v.y() + t.xxz()*v.z(),
        t.xxy()*v.x() + t.xyy()*v.y() + t.xyz()*v.z(),
        t.xxz()*v.x() + t.xyz()*v.y() + t.xzz()*v.z(),
        t.xyy()*v.x() + t.yyy()*v.y() + t.yyz()*v.z(),
        t.xyz()*v.x() + t.yyz()*v.y() + t.yzz()*v.z(),
        t.xzz()*v.x() + t.yzz()*v.y() + t.zzz()*v.z()
    );
}

template<class Cmpt>
class innerProduct<SymmTensor3rdOrder<Cmpt>, Vector<Cmpt>>
{
public:
    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SymmTensor3rdOrder<Cmpt>>
{
public:
    typedef SymmTensor<Cmpt> type;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
