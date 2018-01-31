Simple Artery
=============


Outlet Flow
^^^^^^^^^^^

A simple arterial segment is considered in the first tutorial, subject to a variety of outflow boundary conditions. The input file for the first model is reported here below. ::

  MODEL simpleArtery

  SEGMENT ARTERY 0 20.0 50 0 1 2.0 2.0 14.0 MAT1 NONE 0.0 0 0 FLOW INLETDATA

  DATATABLE INLETDATA LIST
  0.0 14.0 
  1000.0 14.0
  ENDDATATABLE

  SOLVEROPTIONS 0.01 10 1000 2 INLETDATA FLOW 1.0e-6 1 1 TEXT

  MATERIAL MAT1 OLUFSEN 1.06 0.04 1.0 2.0e7 -22.5267 8.65e5

Results
"""""""

The graphs here below illustrate the results:

COMPLETE WITH RESULTS!!!

Outlet Pressure
^^^^^^^^^^^^^^^

In the second example, a constant outlet pressure is applied to the arterial segment. A data table is thus created with a fictious time. The new input file looks like ::

  MODEL simpleArtery

  SEGMENT ARTERY 0 20.0 50 0 1 2.0 2.0 14.0 MAT1 NONE 0.0 0 0 PRESSURE PRESSTABLE

  DATATABLE PRESSTABLE LIST
  0.0 100.0
  ENDDATATABLE

  DATATABLE INLETDATA LIST
  0.0 14.0 
  10.0 14.0
  ENDDATATABLE

  SOLVEROPTIONS 0.01 10 1000 2 INLETDATA FLOW 1.0e-6 1 1 TEXT

  MATERIAL MAT1 OLUFSEN 1.06 0.04 1.0 2.0e7 -22.5267 8.65e5

Results
"""""""

COMPLETE WITH RESULTS!!!

Outlet Resistance
^^^^^^^^^^^^^^^^^

Application of an outlet resistance, results in the following input file ::

  MODEL simpleArtery

  SEGMENT ARTERY 0 20.0 50 0 1 2.0 2.0 0.0 MAT1 NONE 0.0 0 0 RESISTANCE RESTABLE

  DATATABLE RESTABLE LIST
  0.0 100.0
  ENDDATATABLE

  DATATABLE INLETDATA LIST
  0.0 14.0 
  10.0 14.0
  ENDDATATABLE

  SOLVEROPTIONS 0.01 10 1000 4 INLETDATA FLOW 1.0e-6 1 1 TEXT

  MATERIAL MAT1 OLUFSEN 1.06 0.04 1.0 2.0e7 -22.5267 8.65e5

Results
"""""""

COMPLETE WITH RESULTS!!!

Outlet Admittance
^^^^^^^^^^^^^^^^^

Finally, we apply an admittance boundary condition, as follows ::

  MODEL simpleArtery

  SEGMENT ARTERY 0 20.0 50 0 1 2.0 2.0 0.0 MAT1 NONE 0.0 0 0 ADMITTANCE ADMTABLE

  DATATABLE ADMTABLE ADMITTANCE 
  100.0 12.5 0.8 1.1 1.0 0
  ENDDATATABLE

  DATATABLE INLETFLOW LIST
  0.0 2.0
  1000.0 2.0
  ENDDATATABLE

  SOLVEROPTIONS 0.001 10 1000 2 INLETFLOW FLOW 1.0e-6 1 1

  MATERIAL MAT1 OLUFSEN 1.06 0.04 1.0 2.0e7 -22.5267 8.65e5

Results
"""""""

COMPLETE WITH RESULTS!!!

