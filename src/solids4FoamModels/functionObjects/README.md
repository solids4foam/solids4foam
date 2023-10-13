# Function Objects

---

Prepared by Ivan Batistić

---

## Section Aims

- This document describes `solids4foam ` function objects which are not available within the standard `OpenFOAM` package;

- Function objects are various post-processing functionalities that are executing during simulation run time;

- Function objects are placed at the bottom of the `system/controlDict` file:

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
          type    solidForces;
          historyPatch     top;
      }
  }
  ```

---

## `cantileverAnalyticalSolution`

- **Function object purpose**
  To generate the analytical solution fields for the bending slender cantilever problem.
  The analytical solution is taken from [[C.E. Augarde, A.J. Deeks, The use of Timoshenko’s exact solution for a  cantilever beam in adaptive analysis. Finite Elements in Analysis and Design, 44, 2008]](https://www.researchgate.net/publication/30053820_The_use_of_Timoshenko's_exact_solution_for_a_cantilever_beam_in_adaptive_analysis):  
  $$
  \sigma_{xx} = \frac{P(L-x)y}{I}, \qquad \sigma_{yy}=0, \qquad \sigma_{xy}=-\frac{P}{2I}\left( \frac{D^2}{4}-y^2\right),
  $$

$$
u_x = \frac{Py}{6EI} \left((6L-3x)x+(2+\nu)\left(y^2-\frac{D^2}{4}\right)  \right),
$$

$$
u_y = \frac{P}{6EI} \left(3\nu y^2(L-x)+(4+5\nu)\frac{D^2x}{4}+(3L-x)x^2  \right),
$$

​	where  $\nu$ is Poisson's ratio, $E$ is Young modulus, $I$ is the second moment of are of the cross-section, $P$ is applied load and $L$ is length of the beam.
```note
Above analytical solution can be used only if the shearing forces on the ends are distributed according to the same parabolic law as the shearing stress $\tau_{xy}$ and the intensity of the normal forces at the built-in end is proportional to $y$.
```

```warning
The current version of the code assumes a rectangular cross-section with unit width and automatically calculates the second moment of inertia!
```

- __Example of usage__

  ```c++
  functions
  {
      cantileverSolution
      {
          type    cantileverAnalyticalSolution;
  
          E      00e6;
          nu      0.3;
          L 	     10;
          D		0.1;
          P		1e5;
          
          //Optional
          cellDisplacement true;
          pointDisplacement true;
          cellStress true;
          pointStress true;
      }
  }
  ```

- __Arguments__

  -  `P` load applied in the minus $y$ direction at the other end of the beam;
  -  `L`  length of the beam;
  -  `D` depth of the beam;
  -  `E` Young's modulus;
  -  `nu` Poisson's ratio.

- __Optional arguments__

  - `cellDisplacement`   write analytical solution for cell-centred displacement field; default is true;
  - `pointDisplacement`   write analytical solution for vertex-centred displacement field; default is true;
  - `cellStress `    write analytical solution for cell-centred stress field; default is true;
  - `pointStress `  write analytical solution for vertex-centred stress field; default is true;

- __Outputs__

  - Analytical solution for the stress tensor field `analyticalStress` in time directories;
  - Analytical solution for the displacement field `analyticalD` in time directories.
  - `cellStressDifference` field; difference between analytical stress and calculated one: `analyticalStress-sigma`;
  - `DDiference` field; difference between analytical displacement and calculated one: `analyticalD-D`.
  - Log at the end of each time-step:
    `Component: 1`  
    `Norms: mean L1, mean L2, LInfL: 0.12 0.2 0.5 `  
    `...`  

- __Tutorial case in which it is used__
  `solids/linearElasticity/cantilever2d/vertexCentredCantilever2d`

---

## `contactPatchTestAnalyticalSolution`

- **Function object purpose**
  To generate the analytical solution for contact patch test.
  The analytical solution for the stress field is taken from [[ Crisfield MA. Re-visiting the contact patch test. Int J Numer Methods Eng. 2000]](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0207%2820000530%2948%3A3%3C435%3A%3AAID-NME891%3E3.0.CO%3B2-V):  
  $$
  \sigma_{x} = \tau_{xy}=0\qquad \sigma_{y} = \dfrac{E}{1-\nu^2}\Delta \qquad \sigma_z = \nu \sigma_y.
  $$

​	where $E$ is Young's modulus, $\nu$ Poisson's ratio and $\Delta$ prescribed displacement of upper block top surface. 

```note
To use this analytical solution, the bottom surface of the lower block must freely deform in the tangential direction. It can be fixed only in the case of zero Poisson's ratio.
```

- __Example of usage__

  ```c++
  functions
  {
  	analyticalSolution
      {
          type    contactPatchTestAnalyticalSolution;
  
          displacement   0.01;
  
          E       1e6;
          nu      1e-15;
      }
  }
  ```

- __Arguments__

  -  `displacement` upper block prescribed vertical displacement;
  -  `E` Young's modulus;
  -  `nu` Poisson's ratio.

- __Optional arguments__

  - None.

- __Outputs__

  - Analytical solution for stress tensor field `analyticalStress` in time directories.

  - Scalar field of relative error named `relativeError` and defined as:
    $$
    e(\%)=\dfrac{\left| \sigma_y - \sigma_y^{analytical} \right|}{\left|\sigma_y^{analytical}\right|} \cdot 100.
    $$

- __Tutorial case in which it is used__
  `solids/linearElasticity/contactPatchTest`  

---

##  `curvedCantileverAnalyticalSolution`

- **Function object purpose**
  To generate the analytical solution for a curved cantilever beam with end loading and traction-free inner and outer surface.
  The analytical solution for the stress field is taken from [[ Sadd MH. Elasticity: Theory, Applications, and Numerics. Elsevier 2009]](https://www.sciencedirect.com/book/9780123744463/elasticity):  
  $$
  N = {a}^2 - {b}^2 + ({a}^2+{b}^2)\;\text{ln}\left(\frac{b}{a}\right)
  $$

$$
\sigma_{r} = \frac{P}{N}\left(r+\frac{a^2b^2}{r^3}-\frac{a^2+b^2}{r}\right)\sin (\theta),
$$

$$
\sigma_{\theta} = \frac{P}{N}\left(3r-\frac{a^2b^2}{r^3}-\frac{a^2+b^2}{r}\right)\sin (\theta),
$$

$$
\tau_{r\theta} = \tau_{\theta r} =-\frac{P}{N}\left(r+\frac{a^2b^2}{r^3}-\frac{a^2+b^2}{r}\right)\cos (\theta),
$$

​	where $a$ is beam inner radius, $b$ is beam outer radius and $P$ is applied shear force.

- __Example of usage__

  ```c++
  functions
  {
      analyticalSolution
      {
          type    curvedCantileverAnalyticalSolution;
    
          rInner  0.31;
          rOuter  0.33;
   
          force   4;
          E       100;
          nu      0.3;
      }
  }
  ```

- __Arguments__

  -  `rInner` inner beam radius;
  -  `rOuter` outer beam radius;
  -  `force` applied force in (N/m) at beam free end;
  -  `E` Young's modulus;
  -  `nu` Poisson's ratio.

- __Optional arguments__

  - None.

- __Outputs__

  - Analytical solution for stress tensor field `analyticalStress` in time directories.

- __Tutorial case in which it is used__
  `solids/linearElasticity/curvedCantilever`  

----

## `hotCylinderAnalyticalSolution`

- **Function object purpose**
  To generate the analytical solution (temperature and stress) fields for the case of a thermally-stressed pipe/cylinder.
  The analytical solution is taken from [[ Timoshenko, Stephen. *Theory of elasticity*. Oxford, 1951.]](https://asmedigitalcollection.asme.org/appliedmechanics/article/37/3/888/427761/Theory-of-Elasticity-3rd-ed):  
  $$
  \sigma_r = \frac{\alpha E \Delta T}{2(1-\nu)\ln\frac{b}{a}}\left( -\ln \frac{b}{r} - \frac{a^2}{(b^2-a^2)}\left( 1-\frac{b^2}{r^2} \right) \ln \frac{b}{a} \right),
  $$

$$
\sigma_{\theta} = \frac{\alpha E \Delta T}{2(1-\nu)\ln\frac{b}{a}}\left( 1-\ln \frac{b}{r} - \frac{a^2}{(b^2-a^2)}\left( 1+\frac{b^2}{r^2} \right) \ln \frac{b}{a} \right),
$$

$$
T = \displaystyle{\frac{\Delta T}{\ln \frac{b}{a}} \ln \frac{b}{r}},
$$

​	where $a$ is pipe inner radius, $b$ is pipe outer radius,  $\nu$ is Poisson's ratio, $E$ is Young modulus, $\alpha$ is coefficient of linear thermal expansion and $\Delta T$ is 	temperature difference between inner and outer pipe surface.

- __Example of usage__

  ```c++
  functions
  {
      analyticalHotCylinder
      {
          type    hotCylinderAnalyticalSolution;
  
          rInner  0.5;
          rOuter  0.7;
  
          TInner  100;
          TOuter  0;
  
          E       200e9;
          nu      0.3;
          alpha   1e-5;
      }
  }
  ```

- __Arguments__

  -  `rInner` inner pipe radius;
  -  `rOuter` outer pipe radius;
  -  `TInner` temperature on the inner pipe surface;
  -  `TOuter` temperature on the outer pipe surface;
  -  `E` Young's modulus;
  -  `nu` Poisson's ratio;
  -  `alpha` coefficient of linear thermal expansion.

- __Optional arguments__

  - None.

- __Outputs__

  - Analytical solution for hoop stress field `analyticalHoopStress` in time directories.
  - Analytical solution for radial stress field `analyticalRadialStress` in time directories.
  - Analytical solution for temperature field `analyticalT` in time directories.

- __Tutorial case in which it is used__
  `solids/thermoelasticity/hotCylinder/hotCylinder`  

---

## `plateHoleAnalyticalSolution`

- **Function object purpose**
  To generate the analytical solution fields for the "hole in a plate" case.
  The analytical solution for the stress field is taken from [[ Timoshenko, Stephen. *Theory of elasticity*. Oxford, 1951.]](https://asmedigitalcollection.asme.org/appliedmechanics/article/37/3/888/427761/Theory-of-Elasticity-3rd-ed):  
  $$
  \sigma_r = \frac{T}{2}\left( 1-\frac{a^2}{r^2}\right) + \frac{T}{2} \left( 1+\frac{3a^4}{r^4} - \frac{4a^2}{r^2} \right)cos(2\theta),
  $$

$$
\sigma_{\theta} = \frac{T}{2}\left( 1+\frac{a^2}{r^2}\right) - \frac{T}{2} \left( 1+\frac{3a^4}{r^4} \right)\cos(2\theta),
$$

$$
\sigma_{r\theta} =  - \frac{T}{2} \left( 1-\frac{3a^4}{r^4} + \frac{2a^2}{r^2} \right)\sin(2\theta),
$$

​	same in cartesian coordinates:
$$
\sigma_{xx} = T \left( 1-\frac{a^2}{r^2}\left(\frac{3}{2}\cos(2\theta)+\cos(4\theta) \right) + \frac{3}{2}\frac{a^4}{r^4}\cos(4\theta) \right),
$$

$$
\sigma_{yy} = T \left( -\frac{a^2}{r^2}\left(\frac{1}{2}\cos(2\theta)-\cos(4\theta) \right) - \frac{3}{2}\frac{a^4}{r^4}\cos(4\theta) \right),
$$

$$
\sigma_{xy} =  T \left( -\frac{a^2}{r^2}\left(\frac{1}{2}\cos(2\theta)+\sin(4\theta) \right) + \frac{3}{2}\frac{a^4}{r^4}\sin(4\theta) \right).
$$

​	Displacement field in cartesian coordinates:
$$
u_x = \frac{Ta}{8\mu}\left( \frac{r}{a}(\kappa+1)\cos\theta+\frac{2a}{r}\left((1+\kappa)\cos(\theta)+\cos (3\theta)\right)-\frac{2a^3}{r^3}\cos(3\theta)  \right),
$$

$$
u_y = \frac{Ta}{8\mu}\left( \frac{r}{a}(\kappa-3)\sin\theta+\frac{2a}{r}\left((1-\kappa)\sin(\theta)+\sin (3\theta)\right)-\frac{2a^3}{r^3}\sin(3\theta)  \right),
$$

​	where $a$ is hole radius, $T$ is far field traction in $x$ direction, $\nu$ is Poisson's ratio, $\mu$ is shear modulus and $\kappa$ parameter is equal to $3-4\nu$.

- __Example of usage__

  ```c++
  functions
  {
      plateHoleSolution
      {
          type    plateHoleAnalyticalSolution;
    
          holeRadius  1;
          farFieldTractionX  1e6;
   
          E       100;
          nu      0.3;
          
          //Optional
          cellDisplacement true;
          pointDisplacement true;
          cellStress true;
          pointStress true;
      }
  }
  ```

- __Arguments__

  -  `holeRadius` radius of the hole centred on the origin;
  -  `farFieldTractionX`  far-field traction in the $x$ direction;
  -  `E` Young's modulus;
  -  `nu` Poisson's ratio.

- __Optional arguments__

  - `cellDisplacement`   write analytical solution for cell-centred displacement field; default is true;
  - `pointDisplacement`   write analytical solution for vertex-centred displacement field; default is true;
  - `cellStress `    write analytical solution for cell-centred stress field; default is true;
  - `pointStress `  write analytical solution for vertex-centred stress field; default is true;

- __Outputs__

  - Analytical solution for stress tensor field `analyticalStress` in time directories;
  - Analytical solution for the displacement field `analyticalD` in time directories.

- __Tutorial case in which it is used__
  None.

---

## `fsiConvergenceData`

- **Function object purpose**
  Reports the number of outer correctors required at each time-step to reach convergence of the FSI coupling.

- __Example of usage__

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

- __Arguments__

  -  None

- __Optional arguments__

  - `region` name;  the default value is set to `region0`.

- __Outputs__

  - Output file: `postProcessing/0/fsiConvergenceData.dat` ;

  - Output file format:

    ```
    # Time nFsiCorrectors
    0 15
    1 12
    ...
    ```

- __Tutorial case in which it is used__
  `fluidSolidInteraction/heatTransfer/flowOverHeatedPlate`  
  ``fluidSolidInteraction/heatTransfer/thermalCavity`   

---

## `hydrostaticPressure`

- **Function object purpose**
  Outputs the hydrostatic component of the stress tensor field
  $$
  \sigma_h=  -\frac{1}{3}\text{tr}( \mathbf{\sigma}).
  $$
  

  where $\mathbf{\sigma}$ is stress tensor.

- __Example of usage__

  ```c++
  functions
  {
      meanStress
      {
          type    hydrostaticPressure;
      }
  }
  ```

- __Arguments__

  -  None.

- __Optional arguments__

  - None.

- __Outputs__

  - Scalar field `hydrostaticPressure` in time directories;

  - Log at the end of each time-step:

    ```
    Hydrostatic pressure: min = 150, max = 500
    ```

- __Tutorial case in which it is used__
  `None`

---

## `principalStresses`

- **Function object purpose**
  Calculate and write principal stress fields. It assumed that the stress tensor is called `sigma` or `sigmaCauchy`. Three vector fields are created: `sigmaMax`, `sigmaMid`, `sigmaMin`. `sigmaMax` is the most positive/tensile principal stress multiplied by the corresponding principal direction; `sigmaMid` is the middle principal stress multiplied by the corresponding principal direction; `sigmaMin` is the most negative/compressive principal stress multiplied by the corresponding principal direction.

- __Example of usage__

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

- __Arguments__

  -  None.

- __Optional arguments__

  - `region` name;  the default value is set to `region0`.
  - `compressionPositive` specify if compression is considered positive; default is `false`.

- __Outputs__

  -  `sigmaMinDir`  vector field in time directories;
  - `sigmaMin`  scalar field in time directories;
  -  `sigmaMaxDir`  vector field in time directories;
  - `sigmaMax`  scalar field in time directories;
  -  `sigmaMidDir`  vector field in time directories;
  - `sigmaMid`  scalar field in time directories;
  - `sigmaDIff` field; difference between `sigmaMax` and `sigmaMin` fields.

- __Tutorial case in which it is used__
  None

---

## `solidDisplacements`

- **Function object purpose**
  Reports the minimum and maximum values of displacement components together with the arithmetic average value of the displacement.

- __Example of usage__

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

- __Arguments__

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- __Optional arguments__

  - None.

- __Outputs__

  - Output file: `postProcessing/0/solidDisplacements<historyPatch>.dat` ;

  - Output file format:

    ```c++
    # Time minX minY minZ maxX maxY maxZ avX avY avZ
    1 1 1 1 3 4 5 1.2 2 4
    2 1 1 1 3 3 2 1 3 5
    ...
    ```

- __Tutorial case in which it is used__
  `solids/hyperelasticity/longWall`  
  `solids/elastoplasticity/neckingBar`  

## `solidForces`

- **Function object purpose**
  Reports the overall force $\mathbf{f}$  and normal force $f_n$ for specified patch:
  $$
  \mathbf{f} = \sum_{N_f} \mathbf{n}_f \cdot \boldsymbol{\sigma}_f, \qquad {f}_n = \mathbf{f} \cdot \mathbf{n}_f,
  $$
  where $\mathbf{n}_f$ is outward unit normal vector, $\boldsymbol{\sigma}$ is Cauchy stress and $N_f$ is number of faces on specified patch. Subscript $f$ is used to denote face centre value.  In the case of TL formulation, the current boundary unit normal vector $\mathbf{n}_f$ is calculated using total deformation gradient and its Jacobian $J \: \mathbf{F}_{inv}^T \cdot \mathbf{n}_f$.

- __Example of usage__

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

- __Arguments__

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- __Optional arguments__

  - None.

- __Outputs__

  - Output file: `postProcessing/0/solidForces<historyPatch>.dat`;

  - Output file format:

    ```c++
    # Time forceX forceY forceZ normalForce
    1 40 40 50 48
    2 40 60 70 45
    ...
    ```

- __Tutorial case in which it is used__
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
   ``fluidSolidInteraction/3dTube ` 
  `fluidSolidInteraction/3dTubeRobin`
  `fluidSolidInteraction-preCICE/3dTube`

## `solidForcesDisplacements`

- **Function object purpose**
  Reports the overall force $f$ vs arithmetic average displacement $\bar{\mathbf{u}}$ for specified patch:
  $$
  \mathbf{f} = \sum_{N_f} \mathbf{n}_f \cdot \boldsymbol{\sigma}_f,
  $$

  $$
  \bar{\mathbf{u}} = \frac{1}{N_f} \left( \sum_{N_f} \mathbf{u}_f \right),
  $$

  where $\mathbf{n}_f$ is outward unit normal vector, $\boldsymbol{\sigma}$ is Cauchy stress, $\mathbf{u}$ is displacement vector and $N_f$ is number of faces on specified patch. Subscript $f$ is used to denote face centre value.  In the case of TL formulation, the current boundary unit normal vector $\mathbf{n}_f$ is calculated using total deformation gradient and its Jacobian $J \: \mathbf{F}_{inv}^T \cdot \mathbf{n}_f$.

- __Example of usage__

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

- __Arguments__

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- __Optional arguments__

  - None.

- __Outputs__

  - Output file: `postProcessing/0/solidForcesDisplacements<historyPatch>.dat`;

  - Output file format:

    ```c++
    # Time dispX dispY dispZ forceX forceY forceZ
    1 0.1 0.1 0.1 40 40 50
    2 0.1 0.2 0.2 40 60 70
    ...
    ```

- __Tutorial case in which it is used__
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
  E_k = \displaystyle{\frac{1}{2} \sum_{N_P} \rho_P \;( \mathbf{v}_P \cdot \mathbf{v}_P) \; V_p},
  $$
  where $\rho$ is density, $\mathbf{v}$ is velocity, $N_P$ is number of cells and $V$ is volume. Subscript $P$ is used to denote cell centre value.

- __Example of usage__

  ```c++
  functions
  {
      kineticEnergy
      {
          type    solidKineticEnergy;
      }
  }
  ```

- __Arguments__

  -  None.

- __Optional arguments__

  - None.

- __Outputs__

  - Output file: `postProcessing/0/solidKineticEnergy.dat` ;

  - Output file format:

    ```c++
    # Time kineticEnergy
    1 0.10
    2 0.11
    ...
    ```

- __Tutorial case in which it is used__  
  `None`

---

## `solidPointDisplacement`

- **Function object purpose**
  Reports displacement vector value at closest mesh point to the specified point. The closest mesh point is determined using Euclidean distance. The displacement field (defined at cell centres) is interpolated to mesh points using the least squares interpolation.

- __Example of usage__
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

- __Arguments__

  -  `point`  - monitoring point vector.

- __Optional arguments__

  - `region` name;  in the case of structural analysis, the default value is set to `solid` otherwise it is `region0`.

- __Outputs__

  - Output file: `postProcessing/0/solidPointDisplacement_<functionObjectName>.dat` ;

  - Output file format:

    ```c++
    # Time Dx Dy Dz magD
    1 0.10 0.20 0.10 0.24494897
    2 0.11 0.20 0.22 0.26551836
    ...
    ```

- __Tutorial case in which it is used__  
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
  Reports  reports displacement value along a line specified by the user. The displacement field (defined at cell centres) is interpolated to mesh points using the least squares interpolation.

- __Example of usage__

  ```c++
  functions
  {
      pointStress
      {
          type    solidPointDisplacementAlongLine;
          startPoint   (0 0 0);
          endPoint	 (2 0 0);
  
          // Optional
          minDist  1e-5;
          region  "solid";
      }
  }
  ```

```warning
This function object is currently only implemented for serial run!
```

- __Arguments__

  -  `startPoint`  line start point;
  -  `endPoint`  line end point.

- __Optional arguments__

  - `minDist` maximum distance at which mesh point will be included in line plot;
  - `region` name;  in the case of structural analysis, the default value is set to `solid` otherwise it is `region0`.

- __Outputs__

  - Output file: `postProcessing/0/solidPointDisplacementAlongLine<functionObjectName>.dat` ;

  - Output file format:

    ```c++
    # PointID PointCoord Dx Dy Dz mag
    152 (0.2 0.3 0.5) 0.1 0.2 0.1 0.2449
    258 (0.4 0.3 0.5) 0.1 0.2 0.2 0.3
    ...
    ```

- __Tutorial case in which it is used__  
  None.  

## `solidPointStress`

- **Function object purpose**
  Reports stress value at closest mesh point to the specified point. The closest mesh point is determined using Euclidean distance. The stress field (defined at cell centres) is interpolated to mesh points using the least squares interpolation.

- __Example of usage__

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

- __Arguments__

  -  `point`  - monitoring point vector.

- __Optional arguments__

  - `region` name;  in the case of structural analysis, the default value is set to `solid` otherwise it is `region0`.

- __Outputs__

  - Output file: `postProcessing/0/solidPointStress_<functionObjectName>.dat` ;

  - Output file format:

    ```c++
    # Time XX XY XZ YY YZ ZZ
    1 1e6 1e6 1e6 1e6 2e6 3e6 
    2 1e6 3e6 2e6 2e6 2e6 3e6
    ...
    ```

- __Tutorial case in which it is used__  
  `solids/viscoelasticity/viscoTube`  

---

## `solidPointTemperature`

- **Function object purpose**
  Reports temperature value at closest mesh point to the specified point. The closest mesh point is determined using Euclidean distance. The temperature field (defined at cell centres) is interpolated to mesh points using the least squares interpolation.

- __Example of usage__

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

- __Arguments__

  -  `point`  - monitoring point vector.

- __Optional arguments__

  - `region` name;  in the case of structural analysis, the default value is set to `solid` otherwise it is `region0`.

- __Outputs__

  - Output file: `postProcessing/solidPointTemperature_<functionObjectName>.dat` ;

  - Output file format:

    ```c++
    # Time value
    1 225
    2 226
    ...
    ```

- __Tutorial case in which it is used__  
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

  where $\rho$ is density, $\mathbf{u}$ is displacement, $N_P$ is number of cells, $g$ is gravity, $V$ is volume, $\mathbf{r}_{ref}$ is reference point with zero potential energy  and $\mathbf{r}_P$ is positional vector of cell centroid. Subscript $P$ is used to denote cell centre value.

- __Example of usage__

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

- __Arguments__

  - `referencePoint` is a coordinate at which the potential energy is zero.

    ```note
    The value for the uniform gravity field $g$ is specified at `constant/g`!
    ```

- __Optional arguments__

  - None.

- __Outputs__

  - Output file: `postProcessing/0/solidPotentialEnergy.dat` ;

  - Output file format:

    ```c++
    # Time potentialEnergy
    1 500
    2 520
    ...
    ```

- __Tutorial case in which it is used__  
  `None`

```warning
`solidPotentialEnergy` is currently implemented only for linear geometry solid models!
```

---

## `solidStresses`

- **Function object purpose**
  Reports the arithmetic average stress on the patch of a solid:
  $$
  \bar{\boldsymbol{\sigma}} = \frac{1}{N_f} \left( \sum_{N_f} \boldsymbol{\sigma}_f \right)
  $$
  where subscript $f$ is used to denote patch face centre value,  $\boldsymbol{\sigma}$ is Cauchy stress tensor and $N_f$ is overall number of patch faces

- __Example of usage__

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

- __Arguments__

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- __Optional arguments__

  - None.

- __Outputs__

  - Output file: `postProcessing/0/solidStresses<historyPatch>.dat` ;

  - Output file format:

    ```c++
    # Time XX XY XZ YY YZ ZZ
    1 1e6 1e6 1e6 1e6 2e6 3e6 
    2 1e6 3e6 2e6 2e6 2e6 3e6
    ...
    ```

- __Tutorial case in which it is used__  
  `solids/elastoplasticity/cylinderExpansion`
  `solids/hyperelasticity/longWall`

---

## `solidTorque`

- **Function object purpose**
  Reports torque on the specified patch about the given axis:
  $$
  \mathbf{r}_m = (\mathbf{r}_f - \mathbf{r}_{pa}) - \mathbf{a}(\mathbf{a} \cdot (\mathbf{r}_f - \mathbf{r}_{pa})),
  $$

  $$
  \text{Torque} = \sum_{N_f} \mathbf{a} \cdot (\mathbf{r}_m \times (\mathbf{S}_f \cdot \boldsymbol{\sigma})),
  $$

  where $\mathbf{r}_{pa}$ is point on axis, $\mathbf{a}$ is axis direction, $\mathbf{r}_f$ is face centre vector, $\mathbf{S}_f$ boundary face area vector, $\boldsymbol{\sigma}_f$ Cauchy stress vector and $N_f$ number of boundary patch faces. In the case of TL formulation, the current boundary face area vector $\mathbf{S}_f$ is calculated using total deformation gradient and its Jacobian $J \: \mathbf{F}_{inv}^T \cdot \mathbf{S}_f$.

- __Example of usage__

  ```c++
  functions
  {
      patchTorque
      {
          type    solidTorque;
          
          historyPatch     right;
          pointOnAxis      (0 0 0);
          axisDirection    (0 0 1);
          
          //Optional
          stressName	sigma;
      }
  }
  ```

- __Arguments__

  - `pointOnAxis` - point on axis;

  - `axisDirestion` - axis vector, does not require to be normalised;

  - `historyPatch` is the name of the patch.

    ```note
    The non-existing patch name will not stop the simulation.
    ```

- __Optional arguments__

  - `stressName`- stress tensor name,  the default value is set to `sigma`.

- __Outputs__

  - Output file: `postProcessing/0/solidTorque<historyPatch>.dat` ;

  - Output file format:
    ```c++
    # Time torqueX torqueY torqueZ
    1 10 12 13
    2 12 13 13
    ...
    ```

- __Tutorial case in which it is used__
  `None`

---

## `solidTractions`

- **Function object purpose**
  Writes boundary traction as a vector field:
  $$
  \text{Updated Lagrangian (UL):} \qquad \mathbf{t}_b =  \mathbf{n}_b \cdot \boldsymbol{\sigma}
  $$

  $$
  \text{Total Lagrangian (TL):} \qquad \mathbf{t}_b =  \frac{\textbf{F}_{inv}^T \cdot \textbf{n}_b}{|\textbf{F}_{inv}^T \cdot \textbf{n}_b|} \cdot \boldsymbol{\sigma}
  $$

  $$
  \text{Small strain approach:} \qquad \mathbf{t}_b =  \mathbf{n}_b \cdot \boldsymbol{\sigma}
  $$

  where $\boldsymbol{\sigma}$ is Cauchy stress tensor, $\mathbf{n}_b$ boundary outward unit vector and $\mathbf{F}_{inv}$ inverse of total deformation gradient. Note that in case of UL formulation boundary normal $\mathbf{n}_b$ is in current configuration while in the case of TL approach it is in initial configuration.

- __Example of usage__

  ```c++
  functions
  {
      patchTractions
      {
          type    solidTractions;
      }
  }
  ```

- __Arguments__

  -  None.

- __Optional arguments__

  - None.

- __Outputs__

  - Vector field `traction` in time directories;

- __Tutorial case in which it is used__
  `None`

---

## `stressTriaxiality`

- **Function object purpose**
  Outputs the stress triaxiality field (mean i.e. hydrostatic stress divided by equivalent stress):
  $$
  \sigma_h=  -\frac{1}{3}\text{tr}( \mathbf{\boldsymbol{\sigma}}).
  $$

  $$
  T.F. = \frac{-\sigma_h}{\sigma_{eq}}
  $$

  where $T.F$ is triaxiality factor, $\boldsymbol{\sigma}$ is Cauchy stress and $\mathbf{\sigma}_{eq}$ is Von Mises equivalent stress.

- __Example of usage__

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

- __Arguments__

  -  None.

- __Optional arguments__

  - `region` name;  the default value is set to `region0`.

- __Outputs__

  - Scalar field `stressTriaxiality` in time directories;

  - Log at the end of each time-step:

    ```
    Stress triaxiality: min = 0.2, max = 0.4
    ```

- __Tutorial case in which it is used__  
  `None`

---

