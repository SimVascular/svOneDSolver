# Imports
import numpy as np
from mpi4py import MPI
from runParallelSolver import *
from simResultRecord import *

# Enums
ip1D = 0
ip2D = 1
dimType = ip2D

# Set result record
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

# Get input quantities for solver: NEED TO FILL WITH YOUR SPECIFIC SETTINGS!!!
# Folder with the oneDSolver Executable
exePath       = 'insert/exe/folder'
# Folder where the results are to be written
runDIR        = 'insert/test run/folder'
# Location of the input file
inputTemplate = 'insert/template/file'
# Location of the output files containing the result summary
outputFile    = 'insert/result/summary/file'
# Template File for result file e.g., simpleArtery_Flow_ARTERY_ for "simpleArtery" segment
resNameTempl = 'insert/result file/template name'

# Generate input variable realizations
# Generate in first processor and then distribute
if(rank == 0):
  np.random.seed(seed=0)
  totRuns = 15
  if(dimType == ip1D):
    x = np.zeros((15,1))
    x[:,0] = np.random.normal(loc=2.0, scale=0.2**2, size=totRuns)
  elif(dimType == ip2D): 
    x = np.zeros((15,2))
    x[:,0] = np.random.normal(loc=2.0, scale=0.2**2, size=totRuns)
    x[:,1] = np.random.normal(loc=200.0, scale=5.0**2, size=totRuns)
  else:
  	print 'ERROR: Invalid dimtype.'
  	print 'Terminating.'
  	MPI.Finalize()
  	sys.exit(0)
else:
  x = None

# Distribute x to all processors
x = comm.bcast(x,root=0)

# Set record for result interrogation
allRes = []
isRow = 0
rowColNum = 5
res = singleTXTFileExtractor(resNameTempl,ipFlow,isRow,rowColNum)
allRes.append(res)

# Create SolverWrap Object
slv = oneDSolverWrap(exePath,runDIR,inputTemplate,comm,allRes)

# Run Solver
output = slv.run(x)

if(rank == 0):
  # Print Output
  np.savetxt(outputFile,output)
  # Complete
  print 'Solution Completed!'
