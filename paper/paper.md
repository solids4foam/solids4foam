---
title: 'solids4foam-v2.1: A toolbox for performing solid mechanics and fluid-solid interaction simulations in OpenFOAM'
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
2. **Compatibility:** Supports the three main OpenFOAM forks: OpenFOAM.com, OpenFOAM.org, and foam-extend.
3. **Ease of Installation:** The toolbox is easy to install and requires minimal additional dependencies beyond OpenFOAM.


# Features

`solids4foam` employs a modular design, offering generic class interfaces for solid mechanics, fluid dynamics, fluid-solid coupling methods, and solid material models. It also supports all native OpenFOAM modularity, including boundary conditions and function objects.

The `solids4foam-v2.1` release includes the following features:

## Partitioned Fluid-Solid Interaction Coupling Methods

- Fixed under-relaxation [@Tukovic2018a]
- Aitkens accelerated under-relaxation [@Tukovic2018a]
- Interface-quasi-Newton coupling [@Degroote2009]
- Robin-Neumann coupling [@Tukovic2018b]
- Thermo-fluid-solid interaction coupling

## Finite Volume Solid Model Discretizations and Solution Algorithms

- Segregated [@Cardiff2018a], coupled [@Cardiff2016a], and explicit solution algorithms
- Linear geometry (small strain) and nonlinear geometry (finite strain) formulations, including total and updated Lagrangian
- Cell-centered and vertex-centered formulations
- Continuum and plate formulations

## Solid Material Models

- Linear elasticity (isotropic, orthotropic [@Cardiff2014a]), plasticity ($J_2$ [@Cardiff2016b], Mohr-Coulomb [@Tang2015]), viscoelasticity [@Cardiff2018a], thermo-elasticity [@Cardiff2018a], poroelasticity [@Tang2015]
- Hyperelasticity (neo-Hookean, Ogden, Mooney-Rivlin [@Oliveira2022; @Oliveira2023], Fung [@Oliveira2022; @Oliveira2023], Yeoh [@Oliveira2022; @Oliveira2023]), hyperelastoplasticity [@Cardiff2016b]
- Interface to Abaqus material model subroutines (UMATs)

## Solid Boundary Conditions

- Frictional contact (node-to-segment [@Cardiff2012a; @Cardiff2016b], segment-to-segment [@Batistic2022; @Batistic2023])
- Cohesive zone models
- Traction, displacement, rotation

## Fluid Models

- Incompressible (PIMPLE, PIMPLE-overset)
- Multiphase (volume-of-fluid)
- Weakly compressible [@Oliveira2022]

## Function Objects

- Energies, displacements, forces, stresses, principal stresses, torques

## Utilities and Scripts

- Scripts for ensuring compatibility with the main OpenFOAM forks
- Mesh conversion utilities: OpenFOAM to/from Abaqus

## Tutorials

- A suite of example cases and benchmark problems to demonstrate functionality and verify performance

# Software compares to other commonly-used packages



# Acknowledgements

The `solids4foam` toolbox has not directly received funding from any source; nonetheless, many implementations within the toolbox have stemmed from research projects at University College Dublin, the University of Zagreb, and other institutes. With this in mind, Philip Cardiff gratefully acknowledges financial support from the Irish Research Council through the Laureate program, grant number IRCLA/2017/45, and the European Research Council (ERC) under the European Union’s Horizon 2020 Research and Innovation programme (Grant agreement No. 101088740), Science Foundation Ireland (SFI) through I-Form (SFI Grant Numbers 16/RC/3872 and 21/RC/10295 P2), MaREI (SFI Grant Number RC2302 2), NexSys (SFI Grant Number 21/SPP/375), as well as the DJEI/DES/SFI/HEA Irish Centre for High-End Computing (ICHEC) for the provision of computational facilities and support (www.ichec.ie), and the UCD ResearchIT Sonic cluster funded by UCD IT Services and the UCD Research Office.


# References