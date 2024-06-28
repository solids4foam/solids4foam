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
  - name: Ivan Batisti\'{c}
    orcid: 0000-0002-6104-0415
    affiliation: 2
  - name: \v{Z}eljko Tukovi\'{c}
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

`solids4foam` is a toolbox for performing solid mechanics and fluid-solid interaction simulations within the popular OpenFOAM [@OpenFOAM.com; @OpenFOAM.org; @foam-extend] software. The `solids4foam` toolbox offered a rich set of features, including fluid-solid and thermo-fluid-solid interaction coupling algorithms, a suite of solid material models, advanced solid boundary conditions (e.g. frictional contact, fracture), and several discretisations (e.g., cell-centred, vertex-centred) and solution algorithms (coupled, segregated, explicit) for solid mechanics.


# Statement of need

The `solids4foam` toolbox has been designed to address four primary needs:
- There is a desire to solve fluid-solid interactions in OpenFOAM;
- There is a desire to solve advanced solid mechanics problems natively in OpenFOAM;
- There is a need for a modular approach for coupling different solid and fluid procedures in OpenFOAM;
- There is a desire for an extendible framework for research into novel finite volume methods for solid mechanics.

In addition, four principles have been followed in the design of the `solids4foam` toolbox:
- If you can use OpenFOAM, you can use solids4foam;
- The three main OpenFOAM forks are supported by `solids4foam`: OpenFOAM.com, OpenFOAM.org, and foam-extend;
- Easy toolbox should be easy to install and minimise additional dependencies beyond the requirements of OpenFOAM;
- A significant emphasis is placed on code design and code style, following the [OpenFOAM coding style guide](https://openfoam.org/dev/coding-style-guide) closely.




# Acknowledgements

The `solids4foam` toolbox has not directly received funding from any source; nonetheless, many implementations within the toolbox have stemmed from research projects at University College Dublin, the University of Zagreb, and other institutes. With this in mind, Philip Cardiff gratefully acknowledges financial support from the Irish Research Council through the Laureate program, grant number IRCLA/2017/45, and the European Research Council (ERC) under the European Union’s Horizon 2020 Research and Innovation programme (Grant agreement No. 101088740), Science Foundation Ireland (SFI) through I-Form (SFI Grant Number 16/RC/3872), MaREI (SFI Grant Number RC2302 2), NexSys (SFI Grant Number 21/SPP/375), as well as the DJEI/DES/SFI/HEA Irish Centre for High-End Computing (ICHEC) for the provision of computational facilities and support (www.ichec.ie), and the UCD ResearchIT Sonic cluster funded by UCD IT Services and the UCD Research Office.


# References