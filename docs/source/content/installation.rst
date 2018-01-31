Installation
############

Prerequisites
-------------

The following software is required:

- Boost Libraries
- Sphinx  (HTML Manual)
- SWIG Interface generator (Python Interface)
- MPI (mpi4py is used for running parameteric simulations from python)
- SuperLU_MT (if you would like to use this sparse solver)

Building Instructions
---------------------

Generate makefiles with CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Out-of-source is the preferred build modality. To do this, first create a build folder (for example "OneDBin") as a sibling of the main source code folder (probably "OneDSolver") ::

  mkdir OneDBin

Change folder into OneDBin ::

  cd OneDBin

Now you can use the ccmake (command line) or cmake-gui (Qt-based CMake GUI) to create the make files.
For example, using ccmake from the OneDBin folder type ::

  ccmake ../OneDSolver

Configure the code until CMake will let you generate the makefiles. Press "g" (ccmake) to generate the makefiles.

To build the code run (within "OneDBin") ::

  make -j n

where n is the number of processors you want to use to build. 

Generate the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

The sphinx documentation can be generate by typing ::

  make docs

CMake Options
-------------

Three options are available in CMake:

- **buildPy** - Build Python Interface
- **buildDocs** - If ON the documentation is included. it can be built using

::

  make docs

- **sparseSolverType** Use Sparse Matrix Solvers. 

Generate the Python Interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Python libraries can be generated for the oneDSolver by turning on the corresponding cmake option, i.e. ::

  buildPy ON.

A dynamic library will be generated in the **/py** folder. Remeber to include the path to the binary **/py** folder when importing the library. Here's an **example** of calling the oneDSolver python interface ::

  # Imports
  import sys
  sys.path.insert(0, 'path/to/oneDBin/py')
  import numpy as np
  # Import oneDSolver Library
  import oneDSolver as ods

  ods.runOneDSolver("path/to/the/input/file");

This example runs the solver directly from Python, without the need to build and launch the c++ executable.
In the current version a limited number of functionalities is available, i.e., a Python interface is generate for all functions in the main.h file, i.e.
 
- void **WriteHeader** ();
- int  **getDataTableIDFromStringKey** (string key);
- void **createAndRunModel** (cvOneDOptions* opts);
- void **readModelFile** (string inputFile, cvOneDOptions* opts, cvStringVec includedFiles);
- void **readModel** (string inputFile, cvOneDOptions* opts);
- void **runOneDSolver** (string inputFile);

Additional libraries/interfaces will be included in future versions.

Sparse Solver Options
^^^^^^^^^^^^^^^^^^^^^

The discrete linear system of equations can be solved with various methods. 
These has a significant impact on the solver performance and total solution time. 
The following sparse solver options are available:

- **Skyline Solver**. This is a serial skyline matrix solver. 
- **SuperLU_MT Solver**. This is a sparse solver that performs factorization in parallel on a shared memory machine. To use SuperLU_MT you need to download and install its library before building OneDSolver.
- **CSparse Solver**. This is a serial sparse solver. The source codes are included with the OneDSolver source code, so it doesn't require the installation of external libraries. 

Skyline Solver
^^^^^^^^^^^^^^

This is the default solver and does not require any external library to be installed. 
Solutions are typically one order of magniture slower with this solver then other sparse solver options. 

SuperLU_MT Solver
^^^^^^^^^^^^^^^^^

SuperLU_MT is a parallel direct linear solver for shared memory architectures. 
The source code for this library is available at  `this page <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`_.
The build functionality in oneDSolver requires the user to specify the local installation folder for the library using the variable **SUPERLU_DIR**. Let's assume you have already downloaded the library, unzipped its content to, say, "path/to/folder/SuperLU_MT_3.0" and built the library through the make utility. From inside the OneDBin folder (see above) you only have to type the following cmake command ::

  cmake -DsparseSolverType="superlu" -DSUPERLU_DIR="path/to/folder/SuperLU_MT_3.0" ../OneDSolver/

**NOTE**: OneDSolver assumes that your SuperLU_MT library has been built using the **PTHREAD** library. Please follow the instruction below on how to built SuperLU_MT to make sure this is the case.

How to build the SuperLU_MT library so it works with OneDSolver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the documentation in the README file distributed with the SuperLU_MT library. 

Here we assume:

- You are in the "path/to/folder/SuperLU_MT_3.0" folder 
- You are using an Ubuntu operating system or equivalent.

For a different operating system and if you want to link SuperLU_MT to a different multithreading library, please refer to the SuperLU_MT installation documentation.

1. **Delete** your file make.inc.
2. **Copy** the file MAKE_INC/make.pthread one folder up and rename it as make.inc.
3. **Install** the BLAS library on your system and edit the BLASLIB entry in make.inc so that it points to the correct location of your system BLAS library.
4. **Build** the code by typing "make". 

CSparse
^^^^^^^

CSparse is a coincise sparse linear algebra package. 
The source codes are publicly available at `this link <http://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html>`_.
For convenience, these source code have been included in the OneDSolver source code.


