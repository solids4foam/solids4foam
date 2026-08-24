---
sort: 7
---

# Function Objects

---

Prepared by Ivan Batistić with edits by Philip Cardiff

---

## Section Aims

- This document describes `solids4foam` function objects which are not
  available within the standard `OpenFOAM` package;

- Function objects are various post-processing functionalities that are
  executing during simulation run time;

- Function objects are placed at the bottom of the `system/controlDict` file;
  for example:

  ```c++
  functions
  {
      forceDisp
      {
          type          solidForcesDisplacements;
          historyPatch  cylinderFixed;
      }
      patchForce
      {
          type       solidForces;
          historyPatch  top;
      }
  }
  ```

## `fsiConvergenceData`

- **Function object purpose**
  Reports the number of outer correctors required at each time-step to reach
  convergence of the FSI coupling.

- **Example of usage**

  ```c++
  functions
  {
      fsiConvData
      {
          type fsiConvergenceData;
          //Optional
          region fluid;
      }
  }
  ```

- **Arguments**

  - None

- **Optional arguments**

  - `region` name; the default value is set to `region0`.

- **Outputs**

  - Output file: `postProcessing/0/fsiConvergenceData.dat` ;

  - Output file format:

    ```plaintext
    # Time nFsiCorrectors
    0 15
    1 12
    ...
    ```

- **Tutorial case in which it is used:**
  `fluidSolidInteraction/heatTransfer/flowOverHeatedPlate`
  `fluidSolidInteraction/heatTransfer/thermalCavity`

---

## `hydrostaticPressure`

- **Function object purpose**
  Outputs the hydrostatic component of the stress tensor field

  $$
  \sigma_h=  -\frac{1}{3}\text{tr}( \mathbf{\sigma}).
  $$

  where $$\mathbf{\sigma}$$ is stress tensor.

- **Example of usage**

  ```c++
  functions
  {
      meanStress
      {
          type    hydrostaticPressure;
      }
  }
  ```

- **Arguments**

  - None.

- **Optional arguments**

  - None.

- **Outputs**

  - Scalar field `hydrostaticPressure` in time directories;

  - Log at the end of each time-step:

    ```plaintext
    Hydrostatic pressure: min = 150, max = 500
    ```

- **Tutorial case in which it is used:**
  `None`

---

## `principalStresses`

- **Function object purpose**
  Calculate and write principal stress fields. It assumed that the stress tensor
  is called `sigma` or `sigmaCauchy`. Three vector fields are created:
  `sigmaMax`, `sigmaMid`, `sigmaMin`. `sigmaMax` is the most positive/tensile
  principal stress multiplied by the corresponding principal direction;
  `sigmaMid` is the middle principal stress multiplied by the corresponding
  principal direction; `sigmaMin` is the most negative/compressive principal
  stress multiplied by the corresponding principal direction.

- **Example of usage**

  ```c++
  functions
  {
      principalStresses1
      {
          type principalStresses;

          // Optional
          compressionPositive   true;
          region    region0;
      }
  }
  ```

- **Arguments**

  - None.

- **Optional arguments**

  - `region` name; the default value is set to `region0`.
  - `compressionPositive` specify if compression is considered positive; default
    is `false`.

- **Outputs**

  - `sigmaMinDir` vector field in time directories;
  - `sigmaMin` scalar field in time directories;
  - `sigmaMaxDir` vector field in time directories;
  - `sigmaMax` scalar field in time directories;
  - `sigmaMidDir` vector field in time directories;
  - `sigmaMid` scalar field in time directories;
  - `sigmaDIff` field; difference between `sigmaMax` and `sigmaMin` fields.

- **Tutorial case in which it is used:**
  None

---

## `solidDisplacements`

- **Function object purpose**
  Reports the minimum and maximum values of displacement components together
  with the arithmetic average value of the displacement.

- **Example of usage**

  ```c++
  functions
  {
      patchDisplacements
      {
          type    solidDisplacements;
          historyPatch     top;
      }
  }
  ```

- **Arguments**

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidDisplacements<historyPatch>.dat` ;

  - Output file format:

      ```c++
      # Time minX minY minZ maxX maxY maxZ avX avY avZ
      1 1 1 1 3 4 5 1.2 2 4
      2 1 1 1 3 3 2 1 3 5
      ...
      ```

- **Tutorial case in which it is used:**
  `solids/hyperelasticity/longWall`
  `solids/elastoplasticity/neckingBar`

## `solidForces`

- **Function object purpose**
  Reports the overall force $$\mathbf{f}$$ and normal force $$f_n$$ for
  specified patch:

  $$
  \mathbf{f} = \sum_{N_f} \mathbf{n}_f \cdot \boldsymbol{\sigma}_f,
  \qquad {f}_n = \mathbf{f} \cdot \mathbf{n}_f,
  $$

  where $$\mathbf{n}_f$$ is outward unit normal vector, $$\boldsymbol{\sigma}$$
  is Cauchy stress and $$N_f$$ is number of faces on specified patch. Subscript
  $f$ is used to denote face centre value. In the case of TL formulation, the
  current boundary unit normal vector $$\mathbf{n}_f$$ is calculated using total
  deformation gradient and its Jacobian
  $$J \: \mathbf{F}_{inv}^T \cdot \mathbf{n}_f$$.

- **Example of usage**

  ```c++
  functions
  {
      patchForce
      {
          type    solidForces;
          historyPatch     top;
      }
  }
  ```

- **Arguments**

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidForces<historyPatch>.dat`;

  - Output file format:

      ```c++
      # Time forceX forceY forceZ normalForce
      1 40 40 50 48
      2 40 60 70 45
      ...
      ```

- **Tutorial case in which it is used:**
  `solids/linearElasticity/punch`
  `solids/linearElasticity/plateHole`
  `solids/elastoplasticity/pipeCrush`
  `solids/elastoplasticity/uniaxialTension`
  `solids/elastoplasticity/impactBar`
  `solids/elastoplasticity/simpleShear`
  `solids/elastoplasticity/perforatedPlate`
  `solids/elastoplasticity/cylinderCrush`
  `solids/abaqusUMATs/plateHoleTotalDispUMAT`
  `solids/hyperelasticity/plateHoleTotalLag`
  `solids/hyperelasticity/cylinderCrush`
  `fluidSolidInteraction/3dTube`
  `fluidSolidInteraction/3dTubeRobin`
  `fluidSolidInteraction-preCICE/3dTube`

## `solidForcesDisplacements`

- **Function object purpose**
  Reports the overall force $$f$$ vs arithmetic average displacement
  $$\bar{\mathbf{u}}$$ for specified patch:

  $$
  \mathbf{f} = \sum_{N_f} \mathbf{n}_f \cdot \boldsymbol{\sigma}_f,
  $$

  $$
  \bar{\mathbf{u}} = \frac{1}{N_f} \left( \sum_{N_f} \mathbf{u}_f \right),
  $$

  where $$\mathbf{n}_f$$ is outward unit normal vector, $$\boldsymbol{\sigma}$$
  is Cauchy stress, $\mathbf{u}$ is displacement vector and $$N_f$$ is number of
  faces on specified patch. Subscript $$f$$ is used to denote face centre value.
  In the case of TL formulation, the current boundary unit normal vector
  $$\mathbf{n}_f$$ is calculated using total deformation gradient and its
  Jacobian $$J \: \mathbf{F}_{inv}^T \cdot \mathbf{n}_f$$.

- **Example of usage**

  ```c++
  functions
  {
      patchForceDisplacements
      {
          type    solidForcesDisplacements;
          historyPatch     top;
      }
  }
  ```

- **Arguments**

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidForcesDisplacements<historyPatch>.dat`;

  - Output file format:

      ```c++
      # Time dispX dispY dispZ forceX forceY forceZ
      1 0.1 0.1 0.1 40 40 50
      2 0.1 0.2 0.2 40 60 70
      ...
      ```

- **Tutorial case in which it is used:**
  `solids/linearElasticity/punch`
  `solids/linearElasticity/plateHole`
  `solids/elastoplasticity/pipeCrush`
  `solids/elastoplasticity/perforatedPlate`
  `solids/elastoplasticity/cylinderCrush`
  `solids/abaqusUMATs/plateHoleTotalDispUMAT`
  `solids/hyperelasticity/plateHoleTotalLag`

---

## `solidKineticEnergy`

- **Function object purpose**
  Reports the kinetic energy $E_k$ of a solid:

  $$
  E_k = \displaystyle{\frac{1}{2} \sum_{N_P} \rho_P \;( \mathbf{v}_P
  \cdot \mathbf{v}_P) \; V_p},
  $$

  where $$\rho$$ is density, $$\mathbf{v}$$ is velocity, $$N_P$$ is number of
  cells and $$V$$ is volume. Subscript $$P$$ is used to denote cell centre
  value.

- **Example of usage**

  ```c++
  functions
  {
      kineticEnergy
      {
          type    solidKineticEnergy;
      }
  }
  ```

- **Arguments**

  - None.

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidKineticEnergy.dat` ;

  - Output file format:

      ```c++
      # Time kineticEnergy
      1 0.10
      2 0.11
      ...
      ```

- **Tutorial case in which it is used:**
  `None`

---

## `solidPressureMinMax`

- **Function object purpose**
  Reports the minimum and maximum values of the solid hydrostatic pressure:

  $$
  p = -\frac{1}{3}\text{tr}(\boldsymbol{\sigma}).
  $$

  where $$\boldsymbol{\sigma}$$ is the Cauchy stress tensor field `sigma`.

- **Example of usage**

  ```c++
  functions
  {
      pressureMinMax
      {
          type    solidPressureMinMax;
      }
  }
  ```

- **Arguments**

  - None.

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidPressureMinMax.dat` ;

  - Output file format:

      ```c++
      # Time  p_min  p_max  trSigma_min  trSigma_max
      1 150 500 -1500 -450
      2 160 510 -1530 -480
      ...
      ```

- **Tutorial case in which it is used:**
  None

---

## `solidPointDisplacement`

- **Function object purpose**
  Reports displacement vector value at closest mesh point to the specified
  point. The closest mesh point is determined using Euclidean distance. The
  displacement field (defined at cell centres) is interpolated to mesh points
  using the least squares interpolation.

- **Example of usage**

  ```c++
  functions
  {
      pointDisp
      {
          type    solidPointDisplacement;
          point   (105e-6 50e-6 0);
          // Optional
          region  "solid";
      }
  }
  ```

- **Arguments**

  - `point` - monitoring point vector.

- **Optional arguments**

  - `region` name; in the case of structural analysis, the default value is set
    to `solid` otherwise it is `region0`.

- **Outputs**

  - Output file:
    `postProcessing/0/solidPointDisplacement_<functionObjectName>.dat` ;

  - Output file format:

      ```c++
      # Time Dx Dy Dz magD
      1 0.10 0.20 0.10 0.24494897
      2 0.11 0.20 0.22 0.26551836
      ...
      ```

- **Tutorial case in which it is used:**
  `solids/hyperelasticity/cylinderCrush`
  `solids/hyperelasticity/cylindricalPressureVessel`
  `solids/abaqusUMATs/plateHoleTotalDispUMAT`
  `solids/elastoplasticity/perforatedPlate`
  `solids/elastoplasticity/cooksMembrane`
  `solids/viscoelasticity/viscoTube`
  `solids/linearElasticity/cooksMembrane`
  `solids/abaqusUMATs/plateHoleTotalDispUMAT`
  `solids/linearElasticity/wobblyNewton`
  `solids/linearElasticity/plateHole`
  `fluidSolidInteraction/beamInCrossFlow`
  `fluidSolidInteraction/HronTurekFsi3`
  `fluidSolidInteraction/3dTubeRobin`

---

## `solidPointDisplacementAlongLine`

- **Function object purpose**
  Reports reports displacement value along a line specified by the user. The
  displacement field (defined at cell centres) is interpolated to mesh points
  using the least squares interpolation.

- **Example of usage**

  ```c++
  functions
  {
      pointStress
      {
          type    solidPointDisplacementAlongLine;
          startPoint   (0 0 0);
          endPoint  (2 0 0);

          // Optional
          minDist  1e-5;
          region  "solid";
      }
  }
  ```

```warning
This function object is currently only implemented for serial run!
```

- **Arguments**

  - `startPoint` line start point;
  - `endPoint` line end point.

- **Optional arguments**

  - `minDist` maximum distance at which mesh point will be included in line
    plot;
  - `region` name; in the case of structural analysis, the default value is set
    to `solid` otherwise it is `region0`.

- **Outputs**

  - Output file:
    `postProcessing/0/solidPointDisplacementAlongLine<functionObjectName>.dat` ;

  - Output file format:

      ```c++
      # PointID PointCoord Dx Dy Dz mag
      152 (0.2 0.3 0.5) 0.1 0.2 0.1 0.2449
      258 (0.4 0.3 0.5) 0.1 0.2 0.2 0.3
      ...
      ```

- **Tutorial case in which it is used:**
  None.

## `solidPointStress`

- **Function object purpose**
  Reports stress value at closest mesh point to the specified point. The closest
  mesh point is determined using Euclidean distance. The stress field (defined
  at cell centres) is interpolated to mesh points using the least squares
  interpolation.

- **Example of usage**

  ```c++
  functions
  {
      pointStress
      {
          type    solidPointStress;
          point   (0.0075 0 0);
          // Optional
          region  "solid";
      }
  }
  ```

- **Arguments**

  - `point` - monitoring point vector.

- **Optional arguments**

  - `region` name; in the case of structural analysis, the default value is set
    to `solid` otherwise it is `region0`.

- **Outputs**

  - Output file: `postProcessing/0/solidPointStress_<functionObjectName>.dat` ;

  - Output file format:

    ```c++
    # Time XX XY XZ YY YZ ZZ
    1 1e6 1e6 1e6 1e6 2e6 3e6
    2 1e6 3e6 2e6 2e6 2e6 3e6
    ...
    ```

- **Tutorial case in which it is used:**
  `solids/viscoelasticity/viscoTube`

---

## `solidPointTemperature`

- **Function object purpose**
  Reports temperature value at closest mesh point to the specified point. The
  closest mesh point is determined using Euclidean distance. The temperature
  field (defined at cell centres) is interpolated to mesh points using the least
  squares interpolation.

- **Example of usage**

  ```c++
  functions
  {
      pointTemp
      {
          type    solidPointTemperature;
          point   (105e-6 50e-6 0);
          // Optional
          region  "solid";
      }
  }
  ```

- **Arguments**

  - `point` - monitoring point vector.

- **Optional arguments**

  - `region` name; in the case of structural analysis, the default value is set
    to `solid` otherwise it is `region0`.

- **Outputs**

  - Output file: `postProcessing/solidPointTemperature_<functionObjectName>.dat`
    ;

  - Output file format:

      ```c++
      # Time value
      1 225
      2 226
      ...
      ```

- **Tutorial case in which it is used:**
  `None`

---

## `solidPotentialEnergy`

- **Function object purpose**
  Reports the potential energy of a solid:

  $$
  h_P = \frac{\mathbf{g}}{|\mathbf{g}|} \cdot (\mathbf{r}_P+\mathbf{u}_P-\mathbf{r}_{ref}),
  $$

  $$
  E_p = \sum_{N_P} \rho_P \: |\mathbf{g}| \: h_P V_P ,
  $$

  where $$\rho$$ is density, $$\mathbf{u}$$ is displacement, $$N_P$$ is number
  of cells, $$g$$ is gravity, $$V$$ is volume, $$\mathbf{r}_{ref}$$ is reference
  point with zero potential energy and $$\mathbf{r}_P$$ is positional vector of
  cell centroid. Subscript $$P$$ is used to denote cell centre value.

- **Example of usage**

  ```c++
  functions
  {
      potentialEnergy
      {
          type    solidPotentialEnergy;
          referencePoint   (10 50 0);
      }
  }
  ```

- **Arguments**

  - `referencePoint` is a coordinate at which the potential energy is zero.

  ```note
  The value for the uniform gravity field $g$ is specified at `constant/g`!
  ```

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidPotentialEnergy.dat` ;

  - Output file format:

      ```c++
      # Time potentialEnergy
      1 500
      2 520
      ...
      ```

- **Tutorial case in which it is used:**
  `None`

```warning
`solidPotentialEnergy` is currently implemented only for linear geometry solid models!
```

---

## `solidStresses`

- **Function object purpose**
  Reports the arithmetic average stress on the patch of a solid:

  $$
  \bar{\boldsymbol{\sigma}} =
  \frac{1}{N_f} \left( \sum_{N_f} \boldsymbol{\sigma}_f \right)
  $$

  where subscript $$f$$ is used to denote patch face centre value,
  $$\boldsymbol{\sigma}$$ is Cauchy stress tensor and $$N_f$$ is overall number
  of patch faces

- **Example of usage**

  ```c++
  functions
  {
      topStress
      {
         type         solidStresses;
         historyPatch top;
      }
  }
  ```

- **Arguments**

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- **Optional arguments**

  - None.

- **Outputs**

  - Output file: `postProcessing/0/solidStresses<historyPatch>.dat` ;

  - Output file format:

    ```c++
    # Time XX XY XZ YY YZ ZZ
    1 1e6 1e6 1e6 1e6 2e6 3e6
    2 1e6 3e6 2e6 2e6 2e6 3e6
    ...
    ```

- **Tutorial case in which it is used:**
  `solids/elastoplasticity/cylinderExpansion` `solids/hyperelasticity/longWall`

---

## `solidTorque`

- **Function object purpose**
  Reports torque on the specified patch about the given axis:

  $$
  \mathbf{r}_m = (\mathbf{r}_f - \mathbf{r}_{pa}) - \mathbf{a}(\mathbf{a}
  \cdot (\mathbf{r}_f - \mathbf{r}_{pa})),
  $$

  $$
  \text{Torque} = \sum_{N_f} \mathbf{a} \cdot (\mathbf{r}_m
  \times (\mathbf{S}_f \cdot \boldsymbol{\sigma})),
  $$

  where $$\mathbf{r}_{pa}$$ is point on axis, $$\mathbf{a}$$ is axis direction,
  $$\mathbf{r}_f$$ is face centre vector, $$\mathbf{S}_f$$ boundary face area
  vector, $$\boldsymbol{\sigma}_f$$ Cauchy stress vector and $$N_f$$ number of
  boundary patch faces. In the case of TL formulation, the current boundary face
  area vector $$\mathbf{S}_f$$ is calculated using total deformation gradient
  and its Jacobian $$J \: \mathbf{F}_{inv}^T \cdot \mathbf{S}_f$$.

- **Example of usage**

  ```c++
  functions
  {
      patchTorque
      {
          type    solidTorque;

          historyPatch     right;
          pointOnAxis      (0 0 0);
          axisDirection    (0 0 1);
      }
  }
  ```

- **Arguments**

  - `pointOnAxis` - point on axis;

  - `axisDirestion` - axis vector, does not require to be normalised;

  - `historyPatch` is the name of the patch.

      ```note
      The non-existing patch name will not stop the simulation.
      ```

- **Optional arguments**

  - None

- **Outputs**

  - Output file: `postProcessing/0/solidTorque<historyPatch>.dat` ;

- **Tutorial case in which it is used:**
  `None`

---

## `solidTractions`

- **Function object purpose**
  Writes boundary traction as a vector field:

  $$
  \text{Updated Lagrangian (UL):} \qquad \mathbf{t}_b =
  \mathbf{n}_b \cdot \boldsymbol{\sigma}
  $$

  $$
  \text{Total Lagrangian (TL):} \qquad \mathbf{t}_b =
  \frac{\textbf{F}_{inv}^T \cdot \textbf{n}_b}{|\textbf{F}_{inv}^T
  \cdot \textbf{n}_b|} \cdot \boldsymbol{\sigma}
  $$

  $$
  \text{Small strain approach:} \qquad \mathbf{t}_b =  \mathbf{n}_b \cdot \boldsymbol{\sigma}
  $$

  where $$\boldsymbol{\sigma}$$ is Cauchy stress tensor, $$\mathbf{n}_b$$
  boundary outward unit vector and $$\mathbf{F}_{inv}$$ inverse of total
  deformation gradient. Note that in case of UL formulation boundary normal
  $$\mathbf{n}_b$$ is in current configuration while in the case of TL approach
  it is in initial configuration.

- **Example of usage**

  ```c++
  functions
  {
      patchTractions
      {
          type    solidTractions;
      }
  }
  ```

- **Arguments**

  - None.

- **Optional arguments**

  - None.

- **Outputs**

  - Vector field `traction` in time directories;

- **Tutorial case in which it is used:**
  `None`

## `stressTriaxiality`

- **Function object purpose**
  Outputs the stress triaxiality field (mean i.e. hydrostatic stress divided by
  equivalent stress):

  $$
  \sigma_h=  -\frac{1}{3}\text{tr}( \mathbf{\boldsymbol{\sigma}}).
  $$

  $$
  T.F. = \frac{-\sigma_h}{\sigma_{eq}}
  $$

  where $$T.F$$ is triaxiality factor, $$\boldsymbol{\sigma}$$ is Cauchy stress
  and $$\mathbf{\sigma}_{eq}$$ is Von Mises equivalent stress.

- **Example of usage**

  ```c++
  functions
  {
      triaxility
      {
          type    stressTriaxility;
          // Optional
          region  "region0";
      }
  }
  ```

- **Arguments**

  - None.

- **Optional arguments**

  - `region` name; the default value is set to `region0`.

- **Outputs**

  - Scalar field `stressTriaxiality` in time directories;

  - Log at the end of each time-step:

    ```plaintext
    Stress triaxiality: min = 0.2, max = 0.4
    ```

- **Tutorial case in which it is used:**
  `None`

---

## `transformStressToCylindrical`

- **Function object purpose**
  Transform stress tensor to cylindrical coordinate system:

  $$
  \sigma_{transformed} =   \mathbf{R} \cdot \sigma \cdot \mathbf{R}^{T},
  $$

  where $$\mathbf{R}$$ is rotation tensor.

- **Example of usage**

  ```c++
  functions
  {
      transformStressToCylindrical
      {
          type        transformStressToCylindrical;

          origin      (0 0 0);
          axis        (0 0 1);
      }
  }
  ```

- **Arguments**

  - `origin` - origin point;
  - `axis` - axis vector, does not require to be normalised;

- **Optional arguments**

  - None.

- **Outputs**

  - Transformed sigma stress field named `sigma:Transformed` in time
    directories;

    When visualizing in Paraview, keep in mind that the stress components in the
    cylindrical coordinate system will have the same names as the ones from the
    Cartesian coordinate system.

    $$
    \boldsymbol{\sigma}  = \left[\begin{array}{ccc}\sigma_{xx} &
    \sigma_{xy} & \sigma_{xz} \\
    & \sigma_{yy} & \sigma_{yx} \\
    &  & \sigma_{zz}\end{array}\right]\equiv\left[\begin{array}{ccc}\sigma_{R R}
    & \sigma_{R \theta} & \sigma_{R \phi} \\
    & \sigma_{\theta \theta} & \sigma_{\theta \phi} \\
    &  & \sigma_{\phi \phi}\end{array}\right]
    $$

- **Tutorial case in which it is used:**
  `solids/linearElasticity/pressurisedCylinder`
  `solids/thermoelasticity/hotCylinder/hotCylinder`
  `solids/multiMaterial/layeredPipe`

---
