Input File Guide
================

Nomenclature
^^^^^^^^^^^^

In this section of the documentation, the following nomenclature will be used:

* *Card* refers to a command in a line of the input file. 
* A *field* is a single string (no space separators) or number that forms a Card, i.e., a command in the input file.
* A *segment* is a straight segment of vasculature with tapered circular cross section. It connects two joints, one joint and one inlet/outlet or directly one inlet with one outlet. 
* A *joint* is a connection between segments that ensures a unique pressure value and mass continuity.
* A *data table* is a table that can be used to specify scalar and time-varying quantities.

Units
^^^^^

In what follows, we assume that the model quantities and associated boundary conditions are always specified in CGS units. 

MODEL Card
^^^^^^^^^^

This card allows the user to specify a name for the model that will be used when generating the output files. If, for example, we want to call our model *Artery*, we need to specify the following line ::

  MODEL Artery

The following field is required:
1. Model Name (string)

INCLUDE Card
^^^^^^^^^^^^

This card is used to recursively include input files in the project. An example is ::

  INCLUDE auxFile.in TRUE

The following fields are required:

1. Name of the file to include (string, no space separation)
2. Activate file. The following two options are supported: 

  * TRUE. The file is included in the current model
  * FALSE. The file is not used. 


NODE Card
^^^^^^^^^^

This card is used to specify the coordinates of a connection between vessel segments. An example is ::

  NODE 0 1.0 2.0 3.0

The following fields are required:

1. Node Number (double)
2. Node X Coordinate (double)
3. Node Y Coordinate (double)
4. Node Z Coordinate (double)


JOINT Card
^^^^^^^^^^

This card is used to specify a connection between vessel segments. By entering the inlet and outlet vessel segments is it possible to enforce a unique value of pressure in the junction and a flow rate that satisfy conservation of mass. An example is ::

  JOINT JOINT1 1 IN0 OUT0

The following fields are required:

1. Joint Name (string)
2. Joint Node (double)
3. Joint Inlet Name (string)
4. Joint Outlet Name (string)


JOINTINLET Card
^^^^^^^^^^^^^^^

This card is used to specify a list of segments ID numbers as inlets for a joint entity. An example is ::

  JOINTINLET IN0 3 2 4 5

This means that the joint inlet named as *IN0* has 3 inlets with segment ID equal to 2,4 and 5 respectively. 

The following fields are required:

1. Inlet Name (string)
2. Total Number of segments (int)
3. List of segments (list of int)


JOINTOUTLET Card
^^^^^^^^^^^^^^^^

This card is used to specify a list of segments ID numbers as outlets for a joint entity. An example is ::

  JOINTOUTLET OUT0 3 2 4 5

This means that the joint outlet named as *OUT0* has 3 inlets with segment ID equal to 2,4 and 5 respectively. 

The following fields are required:

1. Outlet Name (string)
2. Total Number of segments (int)
3. List of segments (list of int)


SEGMENT Card
^^^^^^^^^^^^

The segment card is used to specify vessel segments. An example is ::

  SEGMENT ARTERY 0 40.0 50 0 1 2.0 2.0 4.0 MAT1 NONE 0.0 0 0 FLOW INLETDATA

The following fields are required:

1. Segment Name (string)
2. Segment ID (int)
3. Segment Length (double)
4. Total Finite Elements in Segment (int)
5. Segment Inlet **Node** (int)
6. Segment Outlet **Node** (int)
7. Segment Inlet Area (double)
8. Segment Outlet Area (double)
9. Segment Initial Flow (double)
10. Segment Material (string)
11. Type of Minor Loss. The following minor pressure losses are supported:

  * *NONE*. No pressure loss. 
  * *STENOSIS*.
  * *BRANCH_THROUGH_DIVIDING*.
  * *BRANCH_SIDE_DIVIDING*.
  * *BRANCH_THROUGH_CONVERGING*.
  * *BRANCH_SIDE_CONVERGING*.
  * *BIFURCATION_BRANCH*.

12. Branch Angle (double)
13. Upstream Segment ID (int)
14. Branch Segment ID (int)
15. Boundary Condition Type. The following boundary conditions are supported:

  * *NOBOUND*. No outlet boundary condition.
  * *PRESSURE*. Constant pressure in the model units. 
  * *AREA*. 
  * *FLOW*. Time-varying outlet flow rate. 
  * *RESISTANCE*. Constant resistance in model units. 
  * *RESISTANCE_TIME*. Time-varying resistance in model units. 
  * *PRESSURE_WAVE*
  * *WAVE*
  * *RCR*. Boundary condition specified through an RCR circuit.
  * *CORONARY*. Coronary boundary condition. 
  * *ADMITTANCE*. Admittance boundary condition.
  * *PULMONARY*. Boundary condition using pulmonary morphometry.

16. Data Table Name for boundary condition (string)

DATATABLE Card
^^^^^^^^^^^^^^

This cards is used to directly specify constant and time-varying quantities for inlet/outlet boundary conditions. It also computes admittance and impedance from a parametric definition of the downstream vessel morphometry. An example with a constant inlet flow rate of 14.0 is::

  DATATABLE INLETDATA LIST
  0.0 14.0 
  10.0 14.0
  ENDDATATABLE

The following fields are required:

1. Data Table Name (string)
2. Data Table Type (string). The following types are supported:

  * *LIST*. List of couples time-values. 
  * *IMPEDANCE*. Computation of impedance 
  * *RCRIMPEDANCE*.
  * *MORPHIMPEDANCE*.
  * *ADMITTANCE*.
  * *RCRADMITTANCE*.
  * *MORPHADMITTANCE*.

3. List of times and Values (e.g., "time0 value0 time1 value1 ..." list of alternating times and values)
4. The card **MUST FINISH** with an ENDDATATABLE command in its own row.

LIST data entries
"""""""""""""""""

If the data table is of type LIST, values are specified by alternating the time and the quantity of interest at that instant in time. An example is ::
  
  DATATABLE TABLE1 LIST
  0.0 0.0 
  1000.0 0.0
  ENDDATATABLE

**Note.** When entering, for example, scalar values for *PRESSURE*, *RESISTANCE*, *RCR*, etc. you need to enter the associated time, even if its value will not be read. The following example shows how to enter an outlet resistance value of 1000.0 Barye s/mL ::

  DATATABLE RTABLE LIST
  0.0 1000.0 
  ENDDATATABLE

IMPEDANCE data entries
""""""""""""""""""""""

The following inputs are required:

1. numTimeSteps = (int)values[0];
2.     lengthRadiusRatio = values[1];
3. Radius at the root of the downstream    rootRadius = values[2];
4.    period = values[3];
5.    scaleFactor = values[4];
6. Number of Fourier modes. Use 0 to obtain values of impedance in time. 

An example is of data table defining an impedance is as follows ::

  DATATABLE IMPDATA IMPEDANCE
  100.0 66.0 0.31 1.1 1.0 0
  ENDDATATABLE

RCRIMPEDANCE data entries
"""""""""""""""""""""""""

The following inputs are required:

1. Number of output steps for impedance.
2. Distal resistance of RCR circuit.
3. Proximal resistance of RCR circuit.
4. Period of cardiac cycle.
5. Capacitance of the RCR circuit.
6. Scale factor that increases the vessel radius for exercise conditions.
7. Number of Fourier modes. Use 0 to obtain values of impedance in time. 

An example is of data table defining an impedance is as follows ::

  DATATABLE IMPDATA IMPEDANCE
  100.0 200.0 10.0 1.1 4.0e-3 1.0 0
  ENDDATATABLE

MORPHIMPEDANCE data entries
"""""""""""""""""""""""""""

The following inputs are required:

1. Number of output steps for impedance.
2. Minimum vessel order in fractal tree.
3. Vessel radius at the vasculature tree root. 
4. Period of the cardiac cycle.
5. Number of Fourier modes. Use 0 to obtain values of impedance in time. 

ADMITTANCE data entries
"""""""""""""""""""""""

The following inputs are required:

1. Number of time steps to solve for.
2. Length to radius ratio of typical downstream vessel.
3. Radius of parent vessel (the fractal tree will be attached to this).
4. Period of cardiac cycle.
5. Scale factor.
6. Integer flag to compute the result in the Fourier (value not zero) or time domain (value equal to zero).

RCRADMITTANCE data entries
""""""""""""""""""""""""""

The following inputs are required:

1. Number of time steps to solve for.
2. Value of RCR distal resistance.
3. Value of RCR proximal resistance.
4. Period of the cardiac cycle.
5. Value of RCR Capacitance.
6. Scale Factor.
7. Integer flag to compute the result in the Fourier (value not zero) or time domain (value equal to zero).

MORPHADMITTANCE data entries
""""""""""""""""""""""""""""

The following inputs are required:

1. Number of output steps for impedance.
2. Minimum vessel order in fractal tree.
3. Vessel radius at the vasculature tree root. 
4. Period of the cardiac cycle.
5. Number of Fourier modes. Use 0 to obtain values of impedance in time otherwise data in the Fourier domain will be returned. 

SOLVEROPTIONS Card
^^^^^^^^^^^^^^^^^^

The SOLVEROPTIONS specifies the option needed by the finite element solver. An example is ::

  SOLVEROPTIONS 0.01 10 1000 4 INLETDATA FLOW 1.0e-3 1 1 TEXT  

The following fields are required:

1. Solver Time Step (double), 
2. Steps Between Saves (int), 
3. Maximum Number of Steps (int)
4. Number of quadrature points for finite elements (int), 
5. Name of Data Table for inlet conditions (string)
6. Type of boundary condition. The following boundary conditions are supported:

  * *NOBOUND*. No outlet boundary condition
  * *PRESSURE*. Constant outlet pressure
  * *AREA*
  * *FLOW*. Time varying outlet flow rate. 
  * *RESISTANCE*. Constant resistance at outlet. 
  * *RESISTANCE_TIME*. Time-varying resistance at the outlet. 
  * *PRESSURE_WAVE*
  * *WAVE*
  * *RCR*. Boundary RCR circuit. 
  * *CORONARY*. Coronary boundary condition. 
  * *ADMITTANCE*. Outlet admittance as boundary condition.  
  * *PULMONARY*. Boundary condition with pulmonary morphometry data. 

7. Convergence tolerance (double)
8. Formulation Type. The following formulations are supported: 

  * *0*. Advective formulation
  * *1*. Conservative formulation

9. Stabilization. The following stabilization options are available:

  * *0*. No stabilization
  * *1*. With stabilization

OUTPUT Card
^^^^^^^^^^^

The OUTPUT card specifies the file formats for the program outputs. An example is ::

  OUTPUT VTK 0

1. Output file format. The following output types are supported:
 
 * TEXT. The output of every segment is written in separate text files for the flow rate, pressure, area and Reynolds number. The rows contain output values at varying locations along the segment while columns contains results at various time instants.
 * VTK. The results for all time steps are plotted to a 3D-like model using the XML VTK file format.
 * BOTH. Both TEXT and VTK results are produced.

2. VTK export option. Two options are available for VTK file outputs (need to specify for BOTH also):
 
 * 0 - Multiple files (default). A separate file is written for each saved increment. A **pvd** file is also provided which contains the time information of the sequence. This is the best option to create animations.
 * 1 - The results for all time steps are plotted to a single XML VTK file.

MATERIAL Card
^^^^^^^^^^^^^

This card is used to specify a constitutive relationship between pressure, cross section diameter and thickness. Example are ::

  MATERIAL MAT1 OLUFSEN 1.06 0.04 1.0 2.0e7 -22.5267 8.65e5
  MATERIAL MAT1 LINEAR 1.06 0.04 1.0 2.0e7

The following fields are required:

1. Material Name (string)
2. Material Type. The following material types are currently supported:

  * *LINEAR*. Linear material model.
  * *OLUFSEN*. Modified material model. 

3. Material Density (double)
4. Material Viscosity (double)
5. Material Reference Pressure (double)
6. Material Exponent (double)
7. Material :math:`k_1` parameter (double)
8. Material :math:`k_2` parameter (double)
9. Material :math:`k_3` parameter (double)

**NOTE**: The reference pressure is the pressure associated with the undeformed area of the vessel. Typically this is the diastolic pressure for a specific vessel. Multiple reference pressures can be assigned to segments using multiple material models.

**NOTE**: For the OLUFSEN material model, all three parameters need to be defined, i.e., :math:`k_1`,:math:`k_2`,:math:`k_3`. For a LINEAR material model instead only the first material parameter :math:`k_1` is used and set equal to :math:`E\,h/r`, i.e., the product of elastic modulus and thickness divided by the radius.
