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
    PetscHandles

Description
    Small Resource Acquisition is Initialization (RAII) helper structs that wrap
    PETSc opaque handles to ensure they are destroyed exactly once and never
    leaked. These wrappers are intentionally minimal: they forbid copying, allow
    moving, and provide `reset()` to destroy the current handle and `put()` to
    obtain an out-parameter suitable for PETSc creation routines, e.g.
    `VecCreateMPI(comm, ..., handles.vec.put())`.

    Currently provided wrappers:
      - VecHandle  : manages a PETSc Vec
      - MatHandle  : manages a PETSc Mat
      - SnesHandle : manages a PETSc SNES
      - KspHandle  : manages a PETSc KSP (useful for custom linear solves)

    Two access functions are provided:
      - get(): returns the raw PETSc pointer (e.g. to a Vec, Mat, SNES or KSP)
      - put(): removes an existing object, clears the pointer, and returns a
               pointer-to-pointer, as required by various PETSc functions

    Example
        VecHandle x;

        // Create a vector
        CHKERRQ( VecCreateMPI(PETSC_COMM_WORLD, nLocal, nGlobal, x.put()) );

        // Use x.get() anywhere a Vec is required
        CHKERRQ( VecSetBlockSize(x.get(), blockSize) );

        // On scope exit (or x.reset()), VecDestroy(&x) is called automatically.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef PetscHandles_H
#define PetscHandles_H

#ifdef USE_PETSC

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>

namespace Foam
{

//- RAII wrapper for PETSc Vec
struct VecHandle
{
    Vec v{nullptr};

    //- default / explicit
    VecHandle() = default;
    explicit VecHandle(Vec in) : v(in) {}

    //- non-copyable
    VecHandle(const VecHandle&) = delete;
    VecHandle& operator=(const VecHandle&) = delete;

    //- movable
    VecHandle(VecHandle&& other) noexcept : v(other.v) { other.v = nullptr; }
    VecHandle& operator=(VecHandle&& other) noexcept
    {
        if (this != &other)
        {
            reset();
            v = other.v;
            other.v = nullptr;
        }
        return *this;
    }

    //- destroy if owned
    ~VecHandle() { reset(); }

    //- destroy and null
    void reset()
    {
        if (v) { VecDestroy(&v); v = nullptr; }
    }

    //- pointer-to-pointer out-param for PETSc create routines
    Vec* put()
    {
        reset();
        return &v;
    }

    //- access
    Vec get() const { return v; }
    explicit operator bool() const { return v != nullptr; }
};

//- RAII wrapper for PETSc Mat
struct MatHandle
{
    Mat m{nullptr};

    MatHandle() = default;
    explicit MatHandle(Mat in) : m(in) {}

    MatHandle(const MatHandle&) = delete;
    MatHandle& operator=(const MatHandle&) = delete;

    MatHandle(MatHandle&& other) noexcept : m(other.m) { other.m = nullptr; }
    MatHandle& operator=(MatHandle&& other) noexcept
    {
        if (this != &other)
        {
            reset();
            m = other.m;
            other.m = nullptr;
        }
        return *this;
    }

    ~MatHandle() { reset(); }

    void reset()
    {
        if (m) { MatDestroy(&m); m = nullptr; }
    }

    Mat* put()
    {
        reset();
        return &m;
    }

    Mat get() const { return m; }
    explicit operator bool() const { return m != nullptr; }
};

//- RAII wrapper for PETSc KSP
struct KspHandle
{
    KSP k{nullptr};

    KspHandle() = default;
    explicit KspHandle(KSP in) : k(in) {}

    KspHandle(const KspHandle&) = delete;
    KspHandle& operator=(const KspHandle&) = delete;

    KspHandle(KspHandle&& other) noexcept : k(other.k) { other.k = nullptr; }
    KspHandle& operator=(KspHandle&& other) noexcept
    {
        if (this != &other)
        {
            reset();
            k = other.k;
            other.k = nullptr;
        }
        return *this;
    }

    ~KspHandle() { reset(); }

    void reset()
    {
        if (k) { KSPDestroy(&k); k = nullptr; }
    }

    KSP* put()
    {
        reset();
        return &k;
    }

    KSP get() const { return k; }
    explicit operator bool() const { return k != nullptr; }
};

//- RAII wrapper for PETSc SNES
struct SnesHandle
{
    SNES s{nullptr};

    SnesHandle() = default;
    explicit SnesHandle(SNES in) : s(in) {}

    SnesHandle(const SnesHandle&) = delete;
    SnesHandle& operator=(const SnesHandle&) = delete;

    SnesHandle(SnesHandle&& other) noexcept : s(other.s) { other.s = nullptr; }
    SnesHandle& operator=(SnesHandle&& other) noexcept
    {
        if (this != &other)
        {
            reset();
            s = other.s;
            other.s = nullptr;
        }
        return *this;
    }

    ~SnesHandle() { reset(); }

    void reset()
    {
        if (s) { SNESDestroy(&s); s = nullptr; }
    }

    SNES* put()
    {
        reset();
        return &s;
    }

    SNES get() const { return s; }
    explicit operator bool() const { return s != nullptr; }
};

} // namespace Foam

#endif // USE_PETSC

#endif // PetscHandles_H

// ************************************************************************* //
