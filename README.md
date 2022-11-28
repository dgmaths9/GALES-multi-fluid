GALES_multi_fluid is an OpenMPI-based parallel finite element numerical code written in C++ to perform multicomponent/multiphase simulations. Simulations can be run on any Linux machine.   



This code repository contains

    1. src directory, which contains the source files
    
    2. sim directory, which has simulation cases in the fluid_mc subdirectory
    
    3. tpl directory, which contains the third-party dependencies to compile and run the simulations
    

Since the code is parallel, OpenMPI must be installed on the machine. Usually, cluster machines are already installed with a suitable version of OpenMPI. Otherwise, it can be built by downloading the source files. On a personal Linux computer, it can be installed through the package repository, for example, on ubuntu by executing "sudo apt install libopenmpi-dev openmpi-bin libhdf5-openmpi-dev". We assume that the OpenMPI is pre-installed on the machine.  

     
The first step is to go into the tpl directory and execute the "build_tpl.sh" file ("build_tpl_ubuntu.sh" if using ubuntu). That will install the code dependencies in the "tpl" directory and add the paths in the ".bashrc" file. 


Compilation: Go to any simulation directory and run the build file (./build) that uses CMake to generate the executable file


Running simulations: The simulation can be run with the mpirun command such as mpirun -np 192 ./executable
