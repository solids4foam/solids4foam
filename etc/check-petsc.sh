#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

die()
{
    echo "PETSc configuration error: $*" >&2
    exit 1
}

find_petsc_header()
{
    local header="$1"
    local candidate
    local candidates=("${PETSC_DIR}/include/${header}")

    if [[ -n "${PETSC_ARCH:-}" ]]; then
        candidates+=("${PETSC_DIR}/${PETSC_ARCH}/include/${header}")
    fi

    for candidate in "${candidates[@]}"; do
        if [[ -f "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done

    return 1
}

read_version_component()
{
    local macro="$1"

    awk -v macro="${macro}" \
        '$1 == "#define" && $2 == macro { print $3; exit }' \
        "${version_header}"
}

if [[ -z "${PETSC_DIR:-}" ]]; then
    echo "NO_PETSC"
    exit 0
fi

[[ -d "${PETSC_DIR}" ]] \
    || die "PETSC_DIR does not name a directory: ${PETSC_DIR}"

config_header="$(find_petsc_header petscconf.h)" \
    || die "could not find petscconf.h under PETSC_DIR=${PETSC_DIR}${PETSC_ARCH:+ with PETSC_ARCH=${PETSC_ARCH}}"
version_header="$(find_petsc_header petscversion.h)" \
    || die "could not find petscversion.h under PETSC_DIR=${PETSC_DIR}${PETSC_ARCH:+ with PETSC_ARCH=${PETSC_ARCH}}"

version_major="$(read_version_component PETSC_VERSION_MAJOR)"
version_minor="$(read_version_component PETSC_VERSION_MINOR)"
version_subminor="$(read_version_component PETSC_VERSION_SUBMINOR)"

if [[ -z "${version_major}" || -z "${version_minor}" || -z "${version_subminor}" ]]; then
    die "could not determine the PETSc version from ${version_header}"
fi

petsc_version="${version_major}.${version_minor}.${version_subminor}"

# Minimum supported PETSc version
#
# PETSc 3.14 and older do not define PETSC_ERR_GPU and PETSC_ERR_MPI, which are
# used in numerics/foamPetscSnesHelper/petscErrorHandling.C, so they cannot
# build solids4foam at all. PETSc 3.15 to 3.24 have been tested with the
# wobblyNewton tutorial using both MUMPS and HYPRE
minimum_version_major=3
minimum_version_minor=15

if (( version_major < minimum_version_major )) \
    || (( version_major == minimum_version_major \
          && version_minor < minimum_version_minor ))
then
    die "PETSc ${petsc_version} is too old: solids4foam requires PETSc ${minimum_version_major}.${minimum_version_minor} or newer"
fi

missing_packages=()

for package in MUMPS HYPRE; do
    if ! grep -Eq \
        "^[[:space:]]*#[[:space:]]*define[[:space:]]+PETSC_HAVE_${package}[[:space:]]+1([[:space:]]|$)" \
        "${config_header}"
    then
        missing_packages+=("${package}")
    fi
done

if (( ${#missing_packages[@]} > 0 )); then
    printf -v missing_list '%s, ' "${missing_packages[@]}"
    missing_list="${missing_list%, }"
    if [[ "${S4F_ALLOW_INCOMPLETE_PETSC:-0}" == "1" ]]; then
        echo "WARNING: PETSc ${petsc_version} is missing ${missing_list}." >&2
        echo "Some solids4foam PETSc tutorials and solver configurations may fail." >&2
    else
        echo "PETSc configuration error: PETSc ${petsc_version} is missing ${missing_list}." >&2
        echo "solids4foam expects PETSc to be configured with MUMPS and HYPRE." >&2
        echo "Reconfigure PETSc with --download-mumps --download-hypre (and their dependencies)," >&2
        echo "or set S4F_ALLOW_INCOMPLETE_PETSC=1 to build without these optional packages." >&2
        exit 1
    fi
else
    echo "Detected PETSc ${petsc_version} with MUMPS and HYPRE." >&2
fi

IFS=' ' read -r config_checksum config_size _ < <(cksum "${config_header}")
IFS=' ' read -r version_checksum version_size _ < <(cksum "${version_header}")

printf 'USE_PETSC|%s|%s|%s|%s:%s|%s:%s\n' \
    "${PETSC_DIR}" \
    "${PETSC_ARCH:-}" \
    "${petsc_version}" \
    "${config_checksum}" \
    "${config_size}" \
    "${version_checksum}" \
    "${version_size}"
