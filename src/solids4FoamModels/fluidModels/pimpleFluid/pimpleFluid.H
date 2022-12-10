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
    pimpleFluid

Description
    pimpleFoam fluid model.

Author
    Based on pimpleFoam.
    Ported by Philip Cardiff, UCD.

SourceFiles
    pimpleFluid.C

\*---------------------------------------------------------------------------*/

// There are significant differences between OpenFOAM/FOAM versions so to
// improve readability and maintainability, the ported versions will be given
// in separate files
#ifdef OPENFOAMFOUNDATION
    #include "pimpleFluid.foundation.H"
#elif OPENFOAMESI
    #include "pimpleFluid.esi.H"
#else
    #include "pimpleFluid.foamextend.H"
#endif

// ************************************************************************* //
