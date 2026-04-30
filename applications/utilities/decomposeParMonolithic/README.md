# `decomposeParMonolithic`

Coupled multi-region decomposition utility for monolithic FSI cases.

## Purpose

`decomposeParMonolithic` partitions multiple volume regions as one graph and
then runs the standard `decomposePar` utility internally with `method manual`
so that a single command produces the normal OpenFOAM processor directory tree
for each region.

This is useful for monolithic FSI cases where one region is much smaller than
the other. Instead of forcing every rank to own both fluid and solid cells, the
combined decomposition can assign some ranks only fluid cells, some only solid
cells, and some both.

The name reflects this target monolithic FSI workflow, but the utility is not
tied to a specific solver class. The generated per-region decomposition can be
used by other multi-region solvers if they support processor ranks with zero
cells in some regions.

## Current workflow

The utility:

1. Reads `system/decomposeParDict` by default, or the file supplied with
   `-decomposeParDict`.
2. Reads the coupled region list and interface definitions from the
   `monolithicCoeffs` sub-dictionary.
3. Loads all listed `fvMesh` regions.
4. Builds one combined graph from:
   - each region's internal cell adjacency
   - extra cross-region edges across the declared FSI interface pairs
5. Runs the selected OpenFOAM decomposition method once on the combined graph.
6. Splits the result back into one `cellDecomposition` list per region.
7. By default, runs `decomposePar` internally for each region using
   `method manual`.

## Dictionary format

Example:

```foam
numberOfSubdomains 4;

method scotch;

monolithicCoeffs
{
    regions (fluid solid);

    regionWeights
    {
        fluid 4;
        solid 3;
    }

    interfaces
    (
        {
            regionA fluid;
            patchA  flag;
            regionB solid;
            patchB  interface;
        }
    );
}
```

### Entries

- `numberOfSubdomains`
  Total MPI ranks for the coupled decomposition.
- `method`
  Standard OpenFOAM decomposition method, for example `scotch`.
- `monolithicCoeffs.regions`
  Volume regions to decompose together.
- `monolithicCoeffs.regionWeights`
  Optional scalar weight per region. The current implementation applies these
  as per-cell weights in the combined graph.
- `monolithicCoeffs.interfaces`
  List of cross-region interfaces used to connect the region graphs.

## Region selection

Region selection is currently dictionary-driven only:

- the coupled region list is read from `monolithicCoeffs.regions`
- command-line `-region` and `-regions` overrides are not currently supported

## Interface assumptions in v1

The current implementation is intentionally narrow:

- it expects conformal interface patches
- it requires the paired patches to have the same number of faces
- it matches the two patches by nearest face centre using an O(n^2) search
- it does not read or interpret a separate mapping mode from the dictionary

In practice, this utility is best suited to direct-map style interfaces such as
the conformal `flag`/`interface` pair in `foilInWind`.

If multiple interface face pairs connect the same two owner cells, the duplicate
cell-cell graph edges are collapsed before decomposition.

## Command-line options

```bash
decomposeParMonolithic
decomposeParMonolithic -decomposeParDict system/decomposeParDict.monolithic
decomposeParMonolithic -cellDist
decomposeParMonolithic -decompose-only
decomposeParMonolithic -force
decomposeParMonolithic -copy-zero
decomposeParMonolithic -no-fields
```

### Supported options

- `-decomposeParDict <file>`
  Use an alternative decompose dictionary.
- `-cellDist`
  Also write `cellDist` fields for visualisation.
- `-decompose-only`
  Only write `cellDecomposition` files and stop before the internal
  `decomposePar` calls.
- `-force`
  Pass `-force` to the internal `decomposePar` calls.
- `-copy-zero`
  Pass `-copyZero` to the internal `decomposePar` calls.
- `-no-fields`
  Pass `-no-fields` to the internal `decomposePar` calls.

## `-decompose-only` mode

`-decompose-only` writes:

- `constant/<region>/cellDecomposition`

for each region and prints the equivalent manual `decomposePar` workflow.

If you complete the decomposition manually afterwards, ensure that your manual
`decomposeParDict` uses the same `numberOfSubdomains` as the monolithic
decomposition that produced the `cellDecomposition` files.

## Output

During the run, the utility prints:

- total cells per region
- matched interface face count per coupled interface
- per-rank cell counts per region
- zero-cell rank counts per region

These statistics are useful for checking whether the coupled decomposition is
doing what you intended.
