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

#include "openFoamTableReader.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "fileOperation.H"
#else
    #include "IFstream.H"
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::openFoamTableReader<Type>::openFoamTableReader(const dictionary& dict)
:
    tableReader<Type>(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::openFoamTableReader<Type>::~openFoamTableReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::openFoamTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, Type>>& data
)
{
    // Read data from file
#ifdef OPENFOAMESIORFOUNDATION
    fileHandler().NewIFstream(fName)()() >> data;
#else
    IFstream(fName)() >> data;
#endif
}


template<class Type>
void Foam::openFoamTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& data
)
{
    // Read data from file
#ifdef OPENFOAMESIORFOUNDATION
    fileHandler().NewIFstream(fName)()() >> data;
#else
    IFstream(fName)() >> data;
#endif
}


// ************************************************************************* //
