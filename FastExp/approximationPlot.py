import matplotlib.pyplot as plt
import numpy as np

fileName = "data.dat"

load    = np.loadtxt(fileName)
x       = load[:,0]
errors  = load[:,1:-1]
expe    = load[:,-1]

plt.figure()
for i in range(len(errors[0,:])) :
	plt.semilogy(x, errors[:,i])
	plt.hold("on")

plt.show()
