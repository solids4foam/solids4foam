/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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
inline SymmTensor4thOrder<Cmpt>::SymmTensor4thOrder()
{}


// Construct given VectorSpace
template<class Cmpt>
inline SymmTensor4thOrder<Cmpt>::SymmTensor4thOrder
(
    const VectorSpace<SymmTensor4thOrder<Cmpt>, Cmpt, 9>& vs
)
:
    VectorSpace<SymmTensor4thOrder<Cmpt>, Cmpt, 9>(vs)
{}


// Construct given three Cmpts
template<class Cmpt>
inline SymmTensor4thOrder<Cmpt>::SymmTensor4thOrder
(
 const Cmpt txxxx, const Cmpt txxyy, const Cmpt txxzz,
                   const Cmpt tyyyy, const Cmpt tyyzz,
                                     const Cmpt tzzzz,
 const Cmpt txyxy,
 const Cmpt tyzyz,
 const Cmpt tzxzx
)
{
    this->v_[XXXX] = txxxx; this->v_[XXYY] = txxyy; this->v_[XXZZ] = txxzz;
                            this->v_[YYYY] = tyyyy; this->v_[YYZZ] = tyyzz;
                                                    this->v_[ZZZZ] = tzzzz;
    this->v_[XYXY] = txyxy;
                            this->v_[YZYZ] = tyzyz;
                                                    this->v_[ZXZX] = tzxzx;
}


// Construct from Istream
template<class Cmpt>
inline SymmTensor4thOrder<Cmpt>::SymmTensor4thOrder(Istream& is)
:
    VectorSpace<SymmTensor4thOrder<Cmpt>, Cmpt, 9>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::xxxx() const
{
    return this->v_[XXXX];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::xxyy() const
{
    return this->v_[XXYY];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::xxzz() const
{
    return this->v_[XXZZ];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::yyyy() const
{
    return this->v_[YYYY];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::yyzz() const
{
    return this->v_[YYZZ];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::zzzz() const
{
    return this->v_[ZZZZ];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::xyxy() const
{
    return this->v_[XYXY];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::yzyz() const
{
    return this->v_[YZYZ];
}

template<class Cmpt>
inline const Cmpt&  SymmTensor4thOrder<Cmpt>::zxzx() const
{
    return this->v_[ZXZX];
}


template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::xxxx()
{
    return this->v_[XXXX];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::xxyy()
{
    return this->v_[XXYY];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::xxzz()
{
    return this->v_[XXZZ];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::yyyy()
{
    return this->v_[YYYY];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::yyzz()
{
    return this->v_[YYZZ];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::zzzz()
{
    return this->v_[ZZZZ];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::xyxy()
{
    return this->v_[XYXY];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::yzyz()
{
    return this->v_[YZYZ];
}

template<class Cmpt>
inline Cmpt& SymmTensor4thOrder<Cmpt>::zxzx()
{
    return this->v_[ZXZX];
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //


//- Double-dot-product between a fourth order symmetric tensor and
//  a symmetric tensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&&(const SymmTensor4thOrder<Cmpt>& t4th, const SymmTensor<Cmpt>& st)
{
  //- s_ij = C_ijkl E_kl
    return SymmTensor<Cmpt>
    (
        t4th.xxxx()*st.xx() + t4th.xxyy()*st.yy() + t4th.xxzz()*st.zz(),
        t4th.xyxy()*st.xy(),
        t4th.zxzx()*st.xz(),

        t4th.xxyy()*st.xx() + t4th.yyyy()*st.yy() + t4th.yyzz()*st.zz(),
        t4th.yzyz()*st.yz(),

        t4th.xxzz()*st.xx() + t4th.yyzz()*st.yy() + t4th.zzzz()*st.zz()
     );
}

//- Double-dot-product between a symmetric tensor and a symmetric
//  fourth order tensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&&(const SymmTensor<Cmpt>& st, const SymmTensor4thOrder<Cmpt>& t4th)
{
  //- s_ij =  E_kl C_klij
  //- as C_ijkl = C_klij
  //- then s_ij = C_ijkl E_kl
    return SymmTensor<Cmpt>
    (
        t4th.xxxx()*st.xx() + t4th.xxyy()*st.yy() + t4th.xxzz()*st.zz(),
        t4th.xyxy()*st.xy(),
        t4th.zxzx()*st.xz(),

        t4th.xxyy()*st.xx() + t4th.yyyy()*st.yy() + t4th.yyzz()*st.zz(),
        t4th.yzyz()*st.yz(),

        t4th.xxzz()*st.xx() + t4th.yyzz()*st.yy() + t4th.zzzz()*st.zz()
    );
}

//- Inner-product between a symmTensor4thOrder and a vector
//  this is not implemented
template<class Cmpt>
inline Vector<Cmpt>
operator&(const SymmTensor4thOrder<Cmpt>& st, const Vector<Cmpt>& v)
{
    notImplemented
      (
       "operator&(const SymmTensor4thOrder<Cmpt>& st, const Vector<Cmpt>& v)"
       " not implemented in SymmTensor4thOrderTemplateI.H"
       );

    return Vector<Cmpt>
    (
        v.x(),
                       v.y(),
                                      v.z()
    );
}


//- Inner-product between a vector and a symmTensor4thOrder
//  this is not implemented
template<class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const SymmTensor4thOrder<Cmpt>& st)
{
    notImplemented
      (
       "operator&(const Vector<Cmpt>& v, const SymmTensor4thOrder<Cmpt>& st)"
       " not implemented in SymmTensor4thOrderTemplateI.H"
       );

    return Vector<Cmpt>
    (
        v.x(),
                       v.y(),
                                      v.z()
    );
}


template<class Cmpt>
class outerProduct<SymmTensor4thOrder<Cmpt>, Cmpt>
{
public:

    typedef SymmTensor4thOrder<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SymmTensor4thOrder<Cmpt> >
{
public:

    typedef SymmTensor4thOrder<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor4thOrder<Cmpt>, Vector<Cmpt> >
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SymmTensor4thOrder<Cmpt> >
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class scalarProduct<SymmTensor4thOrder<Cmpt>, SymmTensor<Cmpt> >
{
public:

  typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class scalarProduct<SymmTensor<Cmpt>, SymmTensor4thOrder<Cmpt> >
{
public:

  typedef SymmTensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
