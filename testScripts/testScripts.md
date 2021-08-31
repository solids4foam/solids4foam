# Test Scripts

## Overview
This directory contains bash scripts to check that solids4foam successfully builds with a number of OpenFOAM/FOAM versions. These scripts are automatically run (by Jenkins) whenever there is a commit pushed to the development branch; if all tests pass then the commit is merged with the master branch.

## How to use the test scripts
**Be careful**: these scripts delete everything from the repository before they finish so it is recommended to only run these scripts on a local copy of the solids4foam repository, which does not contain user cases or code. After running each script, you will need to remove remenants of this copy of solids4foam and re-download it before running another script.

To run a script, download a fresh copy of solids4foam, navigate to the main
solids4foam-release directory and execute the script:
```
$> git clone https://bitbucket.org/philip_cardiff/solids4foam-release.git solids4foam.tmp
$> cd solids4foam.tmp
$> ./testScripts/buildAndTestSolids4foamFE40.sh
$> # Examine the generated *.log files
$> cd ..
$> rm -rf solids4foam.tmp
```

Similarly you can check other versions, e.g.

```
$> git clone https://bitbucket.org/philip_cardiff/solids4foam-release.git solids4foam.tmp
$> cd solids4foam.tmp
$> ./testScripts/buildAndTestSolids4foamOF1912.sh
$> # Examine the generated *.log files
$> cd ..
$> rm -rf solids4foam.tmp
```

**Note**: these scripts require docker to run; they create a temporary docker container to perform the building and testing; this container is then removed at the end of the script. Public docker images are used from Docker Hub, primarily from the philippic account (where you can also find docker images of solids4foam).

**Another Note**: currently scripts check the following:
- Check Allwmake passes without giving an error (FE40, FE41, OF1812, OF1912, OF7)
- Check tutorials/Alltest passes without giving an error (currently only FE40 and FE41)