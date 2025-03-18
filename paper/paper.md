---
title: 'solids4foam: A toolbox for performing solid mechanics and fluid-solid interaction simulations in OpenFOAM'
tags:
  - C++
  - computational mechanics
  - fluid-solid interaction
  - OpenFOAM
  - finite volume method
authors:
  - name: Philip Cardiff
    orcid: 0000-0002-4824-427X
    corresponding: true #
    affiliation: 1
  - name: Ivan Batisti$\'{c}$
    orcid: 0000-0002-6104-0415
    affiliation: 2
  - name: $\v{Z}$eljko Tukovi$\'{c}$
    orcid: 0000-0001-8719-0983
    affiliation: 2
affiliations:
 - name: School of Mechanical and Materials Engineering, University College Dublin, Dublin, Ireland
   index: 1
 - name: Faculty of Mechanical Engineering and Naval Architecture, University of Zagreb, Zagreb, Croatia
   index: 2
date: 28 June 2024
bibliography: paper.bib
---

# Summary

`solids4foam` is a toolbox designed for conducting solid mechanics and fluid-solid interaction simulations within the widely-used OpenFOAM software [@OpenFOAM.com; @OpenFOAM.org; @foam-extend]. The toolbox has a comprehensive set of features, including advanced algorithms for fluid-solid and thermo-fluid-solid coupling, a variety of solid material models, non-trivial solid boundary conditions, and numerous discretisation and solution methods for solid mechanics.


The `solids4foam` toolbox is one of the most comprehensive solid mechanics and fluid-solid interaction toolboxes available within OpenFOAM, having evolved from the `solidMechanics` and `extend-bazaar/FluidSolidInteraction` toolboxes of the foam-extend community [@foam-extend].
Several other OpenFOAM-based toolboxes provide capabilities for solid mechanics and fluid-solid interaction, including `FOAM-FSI` [@Mehl2016], `miniGeotechFoam` [@Tang2015], `explicitSolidDynamics` [@Haider2019, @Haider2017, @Haider2018], as well as more specialised solvers such as the membrane fluid-solid interaction solver [@Wagner2022], a coupled actuator line and finite element analysis tool [@Schmitt2022], and a modular multiphysics framework [@StOnge2023].
However, many of these toolboxes are no longer actively maintained and lack the broad range of solid mechanics and fluid-solid interaction features required for general-purpose simulations. Beyond OpenFOAM-based solutions, preCICE [@Chourdakis2023] provides an alternative approach by coupling OpenFOAM with widely-used finite element solvers such as deal.II [@Arndt2021], CalculiX [@Dhondt2004], FEniCS [@Logg2012], and Code_Aster [@codeAster], enabling flexible multiphysics simulations.
While `solids4foam` is among the first general finite volume-based toolboxes for solid mechanics and fluid-solid interaction, established finite element-based codes such as FEniCS [@Logg2012] and deal.II [@Arndt2021] offer comparable functionality but differ in their numerical methodology.
Furthermore, domain-specific fluid-solid interaction solvers, such as SimVascular [@Zhu2022] and Ambit [@Hirschvogel2024] for cardiovascular simulations or turtlefluid-solid interaction [@Bergersen2020] for general monolithic fluid-solid interaction problems, provide specialised solutions for particular applications.
Despite these alternatives, `solids4foam` remains a uniquely positioned toolbox within OpenFOAM, offering a well-maintained and feature-rich platform for solid mechanics and fluid-solid interaction simulations based on the finite volume method.


# Statement of Need

The `solids4foam` toolbox addresses four primary needs within the OpenFOAM community:

1. The need to perform fluid-solid interactions within OpenFOAM.
2. The need to solve complex solid mechanics problems directly within OpenFOAM.
3. The necessity for a modular approach to coupling various solid and fluid processes in OpenFOAM.
4. The demand for an extendable framework to facilitate research into innovative finite volume methods for solid mechanics.

The design of `solids4foam` adheres to four guiding principles:

1. **Usability:** If you can use OpenFOAM, you can use `solids4foam`.
2. **Compatibility:** Supports the three main OpenFOAM variants: OpenFOAM.com, OpenFOAM.org, and foam-extend.
3. **Ease of Installation:** The toolbox is easy to install and requires minimal additional dependencies beyond OpenFOAM.


# Features

The solids4foam toolbox is designed with a modular architecture, allowing for a flexible and extensible approach to solid mechanics and fluid-solid interaction within OpenFOAM. The core framework consists of generic class interfaces for solid mechanics, fluid dynamics, coupling methods, and solid material models. Additionally, it supports all native OpenFOAM modularity, including boundary conditions and function objects. The `solids4foam-v2.1` release includes the following features:

## Partitioned Fluid-Solid Interaction Coupling Methods

`solids4foam` provides a range of partitioned coupling methods for fluid-solid interaction, including fixed under-relaxation, Aitken’s accelerated under-relaxation, interface-quasi-Newton coupling [@Degroote2009] based on a Dirichlet-Neumann formulation, as well as an added-mass Robin-Neumann formulation. Details of the implementation can be found in [@Cardiff2018a], [@Tukovic2018a], and [@Tukovic2018b]. Thermo-fluid-solid interaction coupling is also available [@solids4foamTFSI].

## Finite Volume Solid Model Discretisations and Solution Algorithms

The toolbox supports multiple discretisation approaches and solution algorithms tailored for solid mechanics. Users can choose between segregated [@Cardiff2018a], coupled [@Cardiff2016a], and explicit solution algorithms, offering flexibility depending on computational constraints and problem requirements. Both linear geometry (small strain) and nonlinear geometry (finite strain) formulations are available, with support for total and updated Lagrangian approaches. `solids4foam` includes both cell-centered and vertex-centered formulations.

## Solid Material Models

A wide range of material models are implemented to cater to various solid mechanics problems. These include linear elasticity, covering isotropic and orthotropic materials [@Cardiff2014a], and inelastic material models such as plasticity (e.g., $J_2$ plasticity [@Cardiff2016b] and Mohr-Coulomb plasticity [@Tang2015]), viscoelasticity [@Cardiff2018a], thermo-elasticity [@Cardiff2018a], and poroelasticity [@Tang2015]. Additionally, the toolbox supports hyperelastic materials, including neo-Hookean, Ogden, Mooney-Rivlin, Fung, and Yeoh models [@Oliveira2022; @Oliveira2023], as well as hyperelastoplasticity [@Cardiff2016b]. To further extend its capabilities, `solids4foam` provides an interface to Abaqus user-defined material subroutines (UMATs).

## Solid Boundary Conditions

To ensure realistic constraints and interactions, solids4foam includes a variety of solid boundary conditions. The toolbox supports frictional contact models, with implementations for node-to-segment [@Cardiff2012a; @Cardiff2016b] and segment-to-segment contact algorithms [@Batistic2022; @Batistic2023]. Additionally, cohesive zone models are available for simulating material failure and debonding processes. Standard boundary conditions for traction, displacement, and rotation constraints are also included.

## Fluid Models

To enable fluid-solid interaction simulations, `solids4foam` integrates with ported versions of the OpenFOAM’s fluid solvers. The toolbox supports incompressible flows, including PIMPLE and PIMPLE-overset methods, as well as multiphase flows using the volume-of-fluid approach. For applications requiring compressibility effects, it also includes a weakly compressible fluid solver [@Oliveira2022], expanding the scope of problems that can be addressed. In addition, `solids4foam` supports coupling via preCICE (see the example in the [tutorials guide](https://www.solids4foam.com/tutorials/)); that is, `solids4foam` can be used as a solid solver in a preCICE-coupled simulation, allowing coupling with a broader range of OpenFOAM and non-OpenFOAM fluid models, potentially with additional physics.

## Function Objects

The toolbox provides several function objects for post-processing and in-situ analysis, allowing users to extract key simulation data. Available function objects include energy calculations, displacement tracking, force evaluation, stress computation, principal stress extraction, and torque measurement.

## Utilities and Scripts

`solids4foam` includes a collection of utilities and scripts to enhance usability and compatibility across different OpenFOAM variants. These utilities ensure smooth interoperability between OpenFOAM versions, simplifying installation and execution. Additionally, solids4foam provides mesh conversion tools to facilitate OpenFOAM-to-Abaqus and Abaqus-to-OpenFOAM mesh transformations, enabling seamless integration with external finite element solvers.

## Tutorials

A set of tutorial cases is included with `solids4foam`, providing users with ready-to-run examples for various solid mechanics and fluid-solid interaction scenarios. These tutorials serve as both educational resources and benchmarks for verifying solver performance, helping users quickly become familiar with the toolbox’s capabilities. As described in the [tutorials guide](https://www.solids4foam.com/tutorials/), the tutorial cases are organised into solids, fluids, fluid-solid interaction, and thermo-fluid-solid interaction cases, with sub-divisions in the solids categories for linear elasticity, elastoplasticity, hyperelasticity, poroelasticity, thermoelasticity, viscoelasticity, multiple materials and fracture.

# Acknowledgements

The `solids4foam` toolbox has not directly received funding from any source; nonetheless, many implementations within the toolbox have stemmed from research projects at University College Dublin, the University of Zagreb, and other institutes. With this in mind, Philip Cardiff gratefully acknowledges financial support from the Irish Research Council through the Laureate program, grant number IRCLA/2017/45, and the European Research Council (ERC) under the European Union’s Horizon 2020 Research and Innovation programme (Grant agreement No. 101088740), Science Foundation Ireland (SFI) through I-Form (SFI Grant Numbers 16/RC/3872 and 21/RC/10295 P2), MaREI (SFI Grant Number RC2302 2), NexSys (SFI Grant Number 21/SPP/375), as well as the DJEI/DES/SFI/HEA Irish Centre for High-End Computing (ICHEC) for the provision of computational facilities and support (www.ichec.ie), and the UCD ResearchIT Sonic cluster funded by UCD IT Services and the UCD Research Office.


# References