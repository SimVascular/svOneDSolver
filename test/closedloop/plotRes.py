import matplotlib.pylab as plt
import numpy as np

# OPTIONS
plotInTime = True

data = np.loadtxt("result.dat")

if(plotInTime):
  size = data.shape[0]
else:
  size = data.shape[1]

for loopA in xrange(0,size):
  if(plotInTime):
    plt.plot(data[loopA,:],label=loopA)
  else:  	
  	plt.plot(data[:,loopA],label=loopA)
  plt.legend()

plt.show()