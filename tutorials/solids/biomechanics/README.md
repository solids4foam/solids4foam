# Biomechanics Tutorials

Solid-mechanics cases for biological tissues using the coupled
pressure-displacement solver (`coupledPressureDisplacementSolid`).

> **Note:** `coupledPressureDisplacementSolid` depends on the `fvBlockMatrix`
> block-coupled linear algebra library, which is currently available only in
> foam-extend. All cases in this directory therefore require foam-extend 4.1.

## cardiac

| Case | Description | Mechanical law |
|---|---|---|
| `heartTissueBeam` | Beam of cardiac tissue under time-dependent pressure | `GuccioneElastic` |
| `ventricleSymm` | Idealised symmetric left ventricle | `GuccioneElastic` |

## vascular

| Case | Description | Mechanical law |
|---|---|---|
| `HGO/ratCarotid` | Rat carotid artery under inflation pressure | `HolzapfelGasserOgdenElastic` |
