# imports 
import sys
import re
import numpy as np

# =====================================
# CLASS TO PROCESS TEMPLATE INPUT FILES
# =====================================
class runFileTemplate():

  # Constructor
  def __init__(self,templateFileName):
    # Read Stream In
    self.stream = open(templateFileName,'r').read()
    parNumbers = re.findall(r'<params,(.*?)>', self.stream)
    parInt = [int(value) for value in parNumbers]
    self.dims = len(np.unique(parInt))

  # Generate Realization File
  def generate(self,x,outFile):
    # Check if the dimensions are compatible
    if(self.dims != len(x)):
      print 'Dimensions in Template: ' + str(self.dims)
      print 'Dimensions in Vector: ' + str(len(x))
      print 'Template file dimensions not compatible. Terminating.'
      sys.exit(-1)
    # Replace Variables
    localStream = self.stream
    for loopA in xrange(len(x)):
      localStream = localStream.replace('<params,' + str(loopA) + '>', '%8.4e' % (x[loopA]))
    # Save to realization file
    text_file = open(outFile,"a+")
    text_file.write(localStream)