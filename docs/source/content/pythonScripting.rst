oneDSolver Python Scripting
###########################

Description of files in the **/py/scripts** folder
--------------------------------------------------

The **/py/scripts** folder contains Python scripts to run **an ensemble of oneDSolver simulations in parallel**.

* **testRunSolver.py** Test routine showing how to run an ensemble of simulations with Python.
* **runParallelSolver.py** File containing the solver wrap classes.
* **simResultRecord.py** Classes to specify the simulation results of interest.
* **fileTemplate.py** Classes to handle file templates.

Running parametric simulations in parallel for optimization and UQ
------------------------------------------------------------------

Preparing a template file
^^^^^^^^^^^^^^^^^^^^^^^^^

For template file preparation, you only need to replace the fixed value of a parameter in the oneDSolver input file with the string **<params,#n>**, where **#n** is the parameter number starting from 0.

An example of a valid template file is ::

  MODEL simpleArtery_Flow_

  SEGMENT ARTERY 0 20.0 50 0 1 <params,0> <params,0> 14.0 MAT1 NONE 0.0 0 0 FLOW INLETDATA

  DATATABLE INLETDATA LIST
  0.0 <params,1>
  1000.0 <params,1>
  ENDDATATABLE

  SOLVEROPTIONS 0.01 10 1000 2 INLETDATA FLOW 1.0e-6 1 1 TEXT

  MATERIAL MAT1 OLUFSEN 1.06 0.04 1.0 2.0e7 -22.5267 8.65e5

The above template contains two parameters. Parameter 0 describes the initial and final segment area, while parameter 1 is used to change the constant inlet flow rate.

**Examples** of template files can be found in the */test* folder with extension **tmpl**.

**Remark**

* When running parametric simulation ensemble, the output mode should be always **TEXT**.

Constructing the wrapper object for oneDSolver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run parametric simulations in Python you need to need to **create an object of the oneDSolverWrap class**, as follows ::

  slv = oneDSolverWrap(exePath,runDIR,inputTemplate,comm,res)

where:

* **exePath** Contains the path to the oneDSolver c++ executable.
* **runDIR** Contains a directory where the simulation files will be stored. A runDIR sub-folder will be created for every ensemble of simulations to be run in parallel. Inside this folder, one directory for each parametric run will contain input files, executable and results.
* **inputTemplate** Specify the location of the template file.
* **comm** Is the MPI communicator (mpi4py is used to access the MPI libraries from Python).
* **res**  Contains an object with the details on how to read the result files (see below).

Specify the results of interest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A **resultRecord** class is used to specify the type of results that should be extracted from each simulation. It is instantiated as :: 

  res = resultRecord(resNameTempl,ipFlow,isRow,rowColNum)

where:

* **resNameTempl** is the result file name template. This name is a combination of the **MODEL** card and **SEGMENT** name. For example for a model definition 

::

  MODEL simpleArtery_Flow_

and a segment named **ARTERY** this quantity should be equal to 

::
  
  resNameTempl = 'simpleArtery_Flow_ARTERY_'

* Result type (integer)

  * **Flow (0)** Volumetric flow rate. 
  * **Area (1)** Segment Area.
  * **Pressure (2)** Blood Pressure. 
  * **Re (3)** Reynolds number of the blood flow.

* **isRow** If this value is set to the integer 0, a **column** of the matrix contained in the result file is extracted. This corresponds to a **single time step and all locations along the segment**. Conversely if an integer 1 is specified, the output will contain the results for a fixed segment location at all time steps.

* **rowColNum** This is the specific number (starting from 0) of the row or column to extract.

Generating numerical parameter realizations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parametric inputs are specified in **matrix form** using a realization matrix. This matrix has:

- Number of rows equal to the number of ensemble simulation you intend to run.

- Number of columns equal to the dimensionality (number of parameter) of your problem. 

**Remarks**

* Make sure **params** is always a matrix. **Do not pass vectors**, even in one-dimension. In this latter case, you should initialize this matrix as (nRun being the number of ensemble simulations needed)

::
  
  params = np.zeros((nRun,1))

* Before running the simulation ensemble you should make sure all the processors have an exact copy of the params vector. You can generate the parameter matrix in the root processor and then use the bcast operation (comm is the MPI communicator in mpi4py)

::
 
  params = comm.bcast(params,root=0) 

Running the simulation ensemble
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once a **oneDSolverWrap** object is instantiated and the parameter realizations generated (e.g., from a d-dimensional quadratue rule in stochastic collocation), running an ensemble of simulations in parallel is easy, you just have to call ::
  
  output = slv.run(params)

The routine returns the outputs **in the same order as inputs in the parameter matrix**. This member function returns a matrix containing

- One set of simulation results in each row. 

- Simulation outputs for a single segment as requiested through the **resultRecord** object.

**Remarks**

- Make sure the number of dimensions in the params matrix **is consistent** with the number of parameters in the template input file.
