# IMPORTS
import os,sys
import math
import numpy as np
from mpi4py import MPI
import uuid
from shutil import copy
import subprocess
from fileTemplate import *

# SET ENUM
ipFlow     = 0
ipArea     = 1
ipPressure = 2
ipRe       = 3

# ====================================================
# SPLIT TOTAL RUNS AMONG PROCESSORS WITH NO EMPTY LOAD
# ====================================================
def splitRuns(size,totRuns):
  procRuns = np.array(np.ones((size,))*math.floor(totRuns/size),dtype=np.int)
  for loopA in xrange(totRuns % size):
    procRuns[loopA] += 1  
  
  # Make a cumulative sum
  procRuns = np.cumsum(procRuns)
  
  # Compute start and end for each processor
  runPtr = np.zeros((size,2),dtype=np.int)
  runPtr[0,0] = 0
  for loopA in xrange(1,size):
    runPtr[loopA,0] = procRuns[loopA - 1]
  for loopA in xrange(size):    
    runPtr[loopA,1] = procRuns[loopA] - 1

  # Return Pointer
  return runPtr

# =========================
# GENERIC SOLVER WRAP CLASS
# =========================
class solverWrap():
  pass

# ===================
# WRAP FOR ONEDSOLVER
# ===================
class oneDSolverWrap(solverWrap):

  # Constructor
  def __init__(self,exePath,runDIR,inputTemplate,comm,resToExtractList):

    # Path to Executable 
    self.exePath = exePath
    # Run folder
    self.runDIR = runDIR
    # Location of Input Template file
    self.inputTemplate = inputTemplate
    # MPI Communicator from mpi4py
    self.comm = comm
    # Information on results to extract
    self.resToExtractList = resToExtractList

  # Gather together results from all processors
  def gatherResults(self,sampleData,numberOfRuns,runPtr):
    maxRunsPerProc = np.max(runPtr[:,1]-runPtr[:,0]+1)
    sendbuf = np.zeros((maxRunsPerProc,),dtype=np.float64)
    # Init Data
    allData = None
    if(self.comm.rank == 0):
      # Assumes all the processors run the same number of simulations
      allData = np.zeros((numberOfRuns,sampleData.shape[1]),dtype=np.float64)
    
    # Loop on the Sample Data Columns
    for loopA in xrange(sampleData.shape[1]):
      
      # Init Data 
      for loopB in xrange(sampleData.shape[0]):
        sendbuf[loopB] = sampleData[loopB,loopA]
      recvbuf = None
      
      # Init Receiving Buffer On Main Processor
      if(self.comm.rank == 0):
        recvbuf = np.zeros((self.comm.size, len(sendbuf)),dtype=np.float64)

      # Gather Data
      self.comm.Gather(sendbuf, recvbuf, root=0)

      # Put together in allData
      if(self.comm.rank == 0):
        # print recvbuf
        for i in xrange(self.comm.size):
          for j in xrange(runPtr[i,0],runPtr[i,1]+1):
            allData[j,loopA] = recvbuf[i,j-runPtr[i,0]]
          
    # Return
    return allData

  # Extract All Results
  def extractResults(self):
    # Extract the result 
    numResQty = len(self.resToExtractList)
    for loopA in xrange(numResQty):
      if(loopA == 0):
        currOut = self.resToExtractList[loopA].getSimulationResults()
        currOut = np.resize(currOut,(len(currOut),1))
      else:
        currOutTmp = self.resToExtractList[loopA].getSimulationResults()
        currOutTmp = np.resize(currOutTmp,(len(currOutTmp),1))
        currOut = np.row_stack((currOut,currOutTmp))    
    # Return Current Outputs
    return np.asarray(currOut)

  # Run OneDSolver in Parallel
  def run(self,x):

    # Not too many processors
    if(self.comm.size > x.shape[0]):
      if(self.comm.rank == 0):
        print 'ERROR: Number of processes exceeds number of runs.'
        print 'Terminating.'
      MPI.Finalize()
      sys.exit(0)

    # Get Total Number of Runs
    numberOfRuns = x.shape[0]

    # Extract Executable Name
    exeName = os.path.basename(self.exePath)

    # The main processor prepares the runs
    if(self.comm.rank == 0):
      print 'Preparation...',
      # Create subdir for this run ensemble
      newpath = self.runDIR + '/' + str(uuid.uuid4())
      if not os.path.exists(newpath):
        os.makedirs(newpath,0777)
      # Create all the folders for the simulations
      for loopA in xrange(numberOfRuns):
        runFLD = newpath + '/' + str(loopA).zfill(4)
        if not os.path.exists(runFLD):
          os.makedirs(runFLD,0777)
        # Copy Executables in folders
        copy(self.exePath,runFLD)
      # Finish Preparation
      print 'OK.'
    else:
      newpath = None
    
    # Distribute newpath string
    newpath = self.comm.bcast(newpath,root=0)
      
    # Wait until the main processor finishes
    self.comm.barrier()

    # Remember the local to global mapping
    runPtr = splitRuns(self.comm.size,x.shape[0])

    # Each Processor takes its part of the realizations
    x = x[runPtr[self.comm.rank,0]:runPtr[self.comm.rank,1]+1,:]

    # Create a template file with parameter realizations and run the simulation
    tmpl = runFileTemplate(self.inputTemplate)
    for loopA in xrange(x.shape[0]):
      # Get Folder Name for this run
      runFLD = newpath + '/' + str(runPtr[self.comm.rank,0] + loopA).zfill(4)

      print 'rank: ',self.comm.rank,' folder: ',runFLD, ' Start: ',runPtr[self.comm.rank,0],' End: ',runPtr[self.comm.rank,1]

      # Generate Template
      tmpl.generate(x[loopA,:],runFLD + '/inputFile.dat')
      # Run Simulations
      os.chdir(runFLD)
      runString = './' + exeName + ' inputFile.dat'     
      shellOut = subprocess.Popen(runString, stdout=subprocess.PIPE, shell=True)
      (out, err) = shellOut.communicate()
      # Get Current Results 
      if(loopA == 0):
        currOutput = self.extractResults()
        currOutput = np.resize(currOutput,(1,len(currOutput)))
      else:
        tmpOutput = self.extractResults()
        tmpOutput = np.resize(tmpOutput,(1,len(tmpOutput)))
        currOutput = np.row_stack((currOutput,tmpOutput))

    # DEBUG: Save single processor file
    # np.savetxt(runFLD + '/outSP_' + str(self.comm.rank),currOutput)
      
    # Wait for all the processors to finish running simulations
    self.comm.barrier()

    # Assemble outputs from different processors
    allData = self.gatherResults(currOutput,numberOfRuns,runPtr)

    # Return Assembled output to main processor
    return allData

