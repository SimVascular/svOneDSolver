import sys,os
import numpy as np

# =============================================
# Get list of files containing a certain string
# =============================================
def getFileContaining(searchString):
  fileList = []
  # Get all the files in the current Dir
  directory = os.listdir('.')
  for fname in directory:
    if searchString in fname:
      fileList.append(fname)
  # Return List
  return fileList

# ======================
# RESULT EXTRACTOR CLASS
# ======================
class resultExtractor():
  pass

# =========================================================
# Class to extract multiple quantities from TXT file result
# =========================================================
class singleTXTFileExtractor(resultExtractor):
  # Constructor
  def __init__(self,fileNameTempl,resultType,isRow,rowColNum):
    # Use tree
    self.fileNameTempl = fileNameTempl
    self.resultType = resultType
    if(isRow == 0):
      self.isRow = False
    else:
      self.isRow = True
    self.rowColNum = rowColNum

  # Construct the result file name
  def getResultFileName(self):
    if(self.resultType == 0):
      # flow
      return self.fileNameTempl + 'flow.dat'
    elif(self.resultType == 1):
      # area
      return self.fileNameTempl + 'area.dat'
    elif(self.resultType == 2):
      # pressure
      return self.fileNameTempl + 'pressure.dat'
    elif(self.resultType == 3):
      # re
      return self.fileNameTempl + 'Re.dat'

  # Get Simulation Results
  def getSimulationResults(self):
    fileName = self.getResultFileName()
    # Open the file with numpy
    resMat = np.loadtxt(fileName)
    if(self.isRow):
      return resMat[self.rowColNum,:]
    else:
      return resMat[:,self.rowColNum]

# =============================================
# Class to extract flows from multiple branches
# =============================================
class multipleFlowExtractor(resultExtractor):
  # Constructor
  def __init__(self,branchNameTmplList):
    # Store local qtys
    self.branchNameTmplList = branchNameTmplList
    self.totBranches = len(branchNameTmplList)
  
  # Get Simulation Results
  def getSimulationResults(self):
    res = np.zeros((self.totBranches,))
    # Loop on the files
    for loopA in xrange(self.totBranches):
      # Get Flow Result File
      fileName = self.branchNameTmplList[loopA] + 'flow.dat'
      # Open the file with numpy
      tmpRead = np.loadtxt(fileName)
      res[loopA] = tmpRead[tmpRead.shape[0]-1,tmpRead.shape[1]-1]
    # Return values
    return np.asarray(res)

# =================================================
# Class to extract Pressures from multiple branches
# =================================================
class multiplePressureExtractor(resultExtractor):
  # Constructor
  def __init__(self,branchNameTmplList):
    # Store local qtys
    self.branchNameTmplList = branchNameTmplList
    self.totBranches = len(branchNameTmplList)
  
  # Get Simulation Results
  def getSimulationResults(self):
    res = np.zeros((self.totBranches,))
    # Loop on the files
    for loopA in xrange(self.totBranches):
      # Get Flow Result File
      fileName = self.branchNameTmplList[loopA] + 'pressure.dat'
      # Open the file with numpy
      tmpRead = np.loadtxt(fileName)
      res[loopA] = tmpRead[tmpRead.shape[0]-1,tmpRead.shape[1]-1]
    # Return values
    return np.asarray(res)

# =======================================================
# Class to extract min and max WSS across the whole model
# =======================================================
class minMaxWSSExtractor(resultExtractor):
  # Constructor
  def __init__(self):
    pass
  
  # Get Simulation Results
  # CAREFULL: Returns the min/max WSS in space at the last timestep!!!
  def getSimulationResults(self):
    maxWSS = -sys.float_info.max
    minWSS = sys.float_info.max
    # Get a list of all files containing WSS
    fileList = getFileContaining('wss.dat')
    # Loop on all the WSS result file
    for loopA in xrange(len(fileList)):
      # Get Flow Result File
      fileName = fileList[loopA]
      # Open the file with numpy
      currBranchWSS = np.loadtxt(fileName)[:,-1]
      maxBranchWSS = currBranchWSS.max()
      minBranchWSS = currBranchWSS.min()
      # Store if needed
      if(maxBranchWSS > maxWSS):
        maxWSS = maxBranchWSS
      if(minBranchWSS < minWSS):
        minWSS = minBranchWSS

    # Return values
    return np.array([minWSS,maxWSS])

    
