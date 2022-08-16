# Installation

The solids4foam toolbox can be installed using one of at least two methods:
- ### **natively**:  with a foam-extend-4.0 installation
- ### **docker**: (www.docker.com)
  - Docker Desktop & Docker Toolbox

Both methods for installing the solids4foam toolbox are explained here.

Automated scripts are provided to test the build as described further in [Test Scripts](../testScripts/testScripts.md).

--------

## [Optional] Dependencies

solids4foam has the same default dependencies as OpenFOAM/foam-extend, so there is no requirement to install anything else. However, if you would like to use the Abaqus UMAT material law functionality then you need to install gfortran. For example, this can be done on Ubuntu with:
```
$> sudo apt install gfortran-5
```
where the version of gfortran should match the version of gcc that is used. In the above example, it is assumed that gcc-5 is used. If using clang or Intel compilers then the version should not matter.

To let solids4foam know that you want to use gfortran, you should set the `S4F_USE_GFORTRAN` environmental variable:
```
$> export S4F_USE_GFORTRAN=1
```

--------

## Installing solids4foam **natively**

- The latest version of solids4foam currently compiles with foam-extend-4.0, foam-extend-4.1, OpenFOAM-v1812, OpenFOAM-v1912, and OpenFOAM-7.

- **Note**: currently all the tutorials only run with the foam-extend-4.0 and foam-extend-4.1 forks. For the other versions, some minor (or not so minor) changes to the tutorials are required: this is currently being addressed!

- One of these forks of OpenFOAM needs to be installed before installing solids4foam, e.g. to install foam-extend-4.0 on your system, following instructions on the OpenFOAM wiki:
https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.0

- ~~To get access to the main development git repository, please email me(philip.cardiff@ucd.ie) with “solids4foam access” in the subject~~. The solids4foam bitbucket repository is now public, and solids4foam can be downloaded using git; it can be placed anywhere convenient:
```
$> cd $FOAM_RUN/..
$> git clone https://bitbucket.org/philip_cardiff/solids4foam-release.git
```

- To build, navigate into the downloaded solids4foam-release directory and run the included Allwmake script:
```
$> cd solids4foam-release
$> ./Allwmake
```

- Once the script finishes, solids4foam should be installed: check by running:
```
$> solids4Foam -help
```

- If everything was installed correctly, you should see the following output:
```
Usage: solids4Foam [-DebugSwitches key1=val1,key2=val2,...] [-DimensionedConstants key1=val1,key2=val2,...] [-InfoSwitches key1=val1,key2=val2,...] [-OptimisationSwitches key1=val1,key2=val2,...] [-Tolerances key1=val1,key2=val2,...] [-case dir] [-dumpControlSwitches] [-noFunctionObjects] [-parallel]  [-help] [-doc] [-srcDoc]
```

solids4Foam is now setup and ready to use!

---

  ## Installing solids4foam using **Docker**

Docker is like a light-weight virtual machine that can allow a Ubuntu (or other) linux installation of solids4foam to run on macOS, Windows or linux. For macOS and Windows, Docker can be installed via **Docker Desktop** or **Docker Toolbox**: Docker Desktop is the preferred approach, while Docker Toolbox is for older Mac and Windows systems that do not meet the requirements of the Docker Desktop.

### **Docker Desktop (Preferred Approach)**

- Install Docker on your system by following the instructions for Windows/macOS/Linux on the Docker website: https://docs.docker.com/install/

  - Select Platform
  - Check System requirements
  - Select "Download from Docker Hub"
  - Select "Get Docker"
  - Complete installation

- Once installed, open a terminal (Powershell terminal on Windows), and pull the solids4foam docker image from Docker Hub (this requires internet connection and the image is ~6.2 GB, it may take a few minutes):
```
$> docker pull philippic/solids4foam-ubuntu18.04
```

- Create a solids4foam container (macOS and linux terminal):
```
$> docker run -itd -v="${HOME}":/home/app/foam/app-4.0/sharedRun --name solids4foam.ubuntu18.04  philippic/solids4foam-ubuntu18.04 /bin/bash
```
- Create a solids4foam container (Windows powershell terminal, where you first need to add your Windows home folder within the docker GUI under Resources -> file sharing -> press the + symbol -> add your home folder):
```
$> docker run -itd -v="${HOME}":/home/app/foam/app-4.0/sharedRun \--name solids4foam.ubuntu18.04 \ philippic/solids4foam-ubuntu18.04 /bin/bash
or (where yourHomeFolderName is replaced with your home folder name)
$> docker run -itd -v=c:/users/yourHomeFolderName:/home/app/foam/app-4.0/sharedRun \--name solids4foam.ubuntu18.04 \ philippic/solids4foam-ubuntu18.04 /bin/bash

```

- You can now attach to the newly created container:
```
$> docker attach solids4foam.ubuntu18.04
$> # you may need to press enter a second time to enter the docker environment
```
- Check the solids4Foam executable exists:
```
$> solids4Foam -help
```
- If everything was installed correctly, you should see the following output:
```
Usage: solids4Foam [-DebugSwitches key1=val1,key2=val2,...] [-DimensionedConstants key1=val1,key2=val2,...] [-InfoSwitches key1=val1,key2=val2,...] [-OptimisationSwitches key1=val1,key2=val2,...] [-Tolerances key1=val1,key2=val2,...] [-case dir] [-dumpControlSwitches] [-noFunctionObjects] [-parallel]  [-help] [-doc] [-srcDoc]
```

- solids4foam is now setup and ready to use! The solids4foam toolbox is installed at /home/app/foam/app-4.0/solids4foam. To exit and stop the docker container, type:
```
$> exit
$> docker stop solids4foam
```
- Note: it is straight-forward to install additional software in the Docker container, e.g.
```
$> sudo apt-get install emacs
```

- Another Note: the shared directory points to the home directory in your host (actual) computer and in the docker container it is located at:
```
/home/app/foam/app-4.0/sharedDir
```

- Once exited, to re-attach to a container that has already been created, you can:
```
$> docker start solids4foam.ubuntu18.04
$> docker attach solids4foam.ubuntu18.04
```

### **Docker Toolbox**

- Install Docker on your system by following the instructions for Windows/macOS/Linux on the Docker website: https://docs.docker.com/toolbox/overview/

  - Select platform for toolbox installation on previous page linked, or follow one of the following links.

    - [Windows](https://docs.docker.com/toolbox/toolbox_install_windows/)
    - [Mac](https://docs.docker.com/toolbox/toolbox_install_mac/)

  - Follow 3 step installation

- Once installed, open a terminal (Docker terminal on Windows), and pull the solids4foam docker image from Docker Hub (this requires internet connection and the image is ~6.2 GB):
```
$> docker pull philippic/solids4foam-ubuntu18.04
```

- Create a solids4foam container (macOS and linux terminal):
```
$> docker run -itd -v="${HOME}":/home/app/foam/app-4.0/sharedRun --name solids4foam.ubuntu18.04  philippic/solids4foam-ubuntu18.04 /bin/bash
```
- Create a solids4foam container (Windows powershell terminal):
```
$> docker run -itd -v="${HOME}":/home/app/foam/app-4.0/sharedRun \--name solids4foam.ubuntu18.04 \ philippic/solids4foam-ubuntu18.04 /bin/bash
```

- You can now attach to the newly created container:
```
$> docker attach solids4foam.ubuntu18.04
$> # you may need to press enter a second time to enter the docker environment
```
- Check the solids4Foam executable exists:
```
$> solids4Foam -help
```
- If everything was installed correctly, you should see the following output:
```
Usage: solids4Foam [-DebugSwitches key1=val1,key2=val2,...] [-DimensionedConstants key1=val1,key2=val2,...] [-InfoSwitches key1=val1,key2=val2,...] [-OptimisationSwitches key1=val1,key2=val2,...] [-Tolerances key1=val1,key2=val2,...] [-case dir] [-dumpControlSwitches] [-noFunctionObjects] [-parallel]  [-help] [-doc] [-srcDoc]
```

- solids4foam is now setup and ready to use! The solids4foam toolbox is installed at /home/app/foam/app-4.0/solids4foam. To exit and stop the docker container, type:
```
$> exit
$> docker stop solids4foam
```
- Note: it is straight-forward to install additional software in the Docker container, e.g.
```
$> sudo apt-get install emacs
```

- Note2: the shared directory points to the home directory in your host (actual) computer and in the docker container it is located at:
```
/home/app/foam/app-4.0/sharedDir
```

- Once exited, to re-attach to a container that has already been created, you can:
```
$> docker start solids4foam.ubuntu18.04
$> docker attach solids4foam.ubuntu18.04
```
---
