### solids4foam ###
solids4foam - a finite volume toolbox for solid mechanics and fluid solid
interaction simulations


### What is this repository for? ###

solids4foam is toolbox for OpenFOAM with capabilities for solid mechanics and
fluid solid interactions.


### How do I get set up? ###

To get setup, you must first install foam-extend-4.0 or foam-extend-4.1,
details of which can be found on the OpenFOAM Wiki at:
https://openfoamwiki.net/index.php/Installation

Once foam-extend-4.0 or foam-extend-4.1 has been installed, download solids4foam
and then run the enclosed `Allwmake` script to compile solids4foam.

solids4foam also compiles with OpenFOAM-7 and OpenFOAM-v1812 but some features
are yet to be ported (porting of all major features is actively underway, Feb-20)
Note that all tutorials have not yet been configured for OpenFOAM-7 and
OpenFOAM-v1812.

Notes:

  1. the `master` branch compiles with foam-extend-4.0 and foam-extend-4.1.
  2. the `master` branch also compiles with OpenFOAM-7 and OpenFOAM-v1812.
  3. the `foam-extend-3.2` branch compiles with foam-extend-3.2 but is not
up-to-date.


### Contribution guidelines ###

If you would like to contribute code and/or test cases to solids4foam, please
email: philip.cardiff@ucd.ie


### Who do I talk to? ###

solids4foam is developed by Philip Cardiff and Zeljko Tukovic, with
contributions for many others, in particular Danial Khazaei.

Emails: philip.cardiff@ucd.ie and zeljko.tukovic@fsb.hr


### Where can I find more information? ###

A number of the tutorial cases are described in the following publications:

P. Cardiff, A Karac, P. De Jaeger, H. Jasak, J. Nagy, A. Ivanković, Ž. Tuković:
An open-source finite volume toolbox for solid mechanics and fluid-solid
interaction simulations. arXiv:1808.10736v2, 2018, available at
https://arxiv.org/abs/1808.10736.

Ž. Tuković, A. Karač, P. Cardiff, H. Jasak, A. Ivanković: OpenFOAM finite volume
solver for fluid-solid interaction.  Transactions of Famena, 42 (3), pp. 1-31,
2018, 10.21278/TOF.42301.

P. Cardiff, A Karac, P. De Jaeger, H. Jasak, J. Nagy, A. Ivanković, Ž. Tuković:
Towards the Development of an Extendable Solid Mechanics and Fluid-Solid
Interactions Toolbox for OpenFOAM. 12th OpenFOAM Workshop University of Exeter,
Exeter, UK. 24th to 27th July 2017.

P. Cardiff, Ž. Tuković, H. Jasak, A. Ivanković: A block-coupled finite vol
me methodology for linear elasticity and unstructured meshes. Computers and
Structures, 2016, 175 100-122, DOI: 10.1016/- j.compstruc.2016.07.004.

P. Cardiff, Ž. Tuković, P. De Jaeger, M. Clancy, A. Ivanković: A Lagrang
an cell-centred finite volume method for metal forming simulation. International
Journal for Numerical Methods in Engineering, 2016, DOI: 10.1002/nme.5345.

P. Cardiff, A. Karać, A. Ivanković: A Large Strain Finite Volume Method f
r Orthotropic Bodies with General Material Orientations. Computer Methods in
Applied Mechanics and Engineering, 01/2014, 268(1):318-335.
DOI: 10.1016/j.cma.2013.09.008.

P. Cardiff, A. Karać, A. Ivanković: Development of a Finite Volume contact solv
r based on the penalty method. Computational Materials Science, 03/2014,
64:283–284. DOI: 10.1016/j.commatsci.2012.03.011.

T. Tang, O. Hededal, P. Cardiff, On Finite Volume method implementation of poro-
elasto-plasticity soil model. International Journal for Numerical and Analytical
Methods in Geomechanics, 2015, DOI: 10.1002/nag.2361.

Ž. Tuković, P. Cardiff, A. Karač, H. Jasak, A. Ivanković: OpenFOAM Library
for Fluid-Structure Interaction. 9th OpenFOAM Workshop, University of Zagreb,
Croatia, 06/2014.
