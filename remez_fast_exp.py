import matplotlib.pyplot as plt
import mpmath as mpm	
import numpy as np
import sys

# Allow the polynomial degree to be given on the command line.
if len(sys.argv) > 1 :
	polynomialDegree = int(sys.argv[1])
else :
	polynomialDegree = 9
N = polynomialDegree

# Function to find minimax polynomial fit for, and its two first derivatives.
f   = lambda x : 1 + x - 2**x
df  = lambda x :  	 1 - 2**x * np.log(2)
ddf = lambda x : 		 2**x * np.log(2)**2

# Start out with a Chebyshev approximation.
poly = mpm.chebyfit(f, [0,1], polynomialDegree+1)
poly = [float(p) for p in poly]

x = np.linspace(0, 1, 5000)
chebyError = np.polyval(poly, x) - f(x)
#plt.plot(x, chebyError, 'r-')
plt.hold('on')

# Find the zeros of the Chebyshev error.
dpoly 	= np.polyder(poly)
dpdf    = lambda x : np.polyval(dpoly, x) - df(x)
xPoints = np.where(np.diff(np.sign(dpdf(x))))[0]
xPoints = x[xPoints]
xPoints = np.insert(xPoints,   0, 0)
xPoints = np.insert(xPoints, N+1, 1)

minE = min(chebyError)
maxE = max(chebyError)
for i in range(polynomialDegree+2) :
	pass
	#plt.plot([xPoints[i], xPoints[i]], [minE, maxE], 'b-')

# Go on to find the minimax approximating polynomial using Remez's algorithm.
def RemezAlgorithm(function, xPoints, degree) :
	A 	= np.zeros((degree+2, degree+2))
	b 	= function(xPoints)
	for i in range(len(A[:,0])) :
		for j in range(len(A[0,:])-1) :
			A[i,j] = xPoints[i]**j
		A[i,-1] = (-1)**(i)
	P 	= np.linalg.solve(A,b)
	eps = P[-1]
	P 	= P[:-1]
	P 	= P[::-1]
	return P, eps

P 				= []
eps 			= 1
maxIterations 	= 1
for iteration in range(maxIterations) :
	if eps < 1e-8 :
		print("HEI")
		break
	P, eps = RemezAlgorithm(f, xPoints, polynomialDegree)

	dP  = np.polyder(P)
	ddP = np.polyder(dP)

	dG  = lambda x : np.polyval(dP, x)  -  df(x)
	ddG = lambda x : np.polyval(ddP, x) - ddf(x)

	for j in range(len(xPoints)) :
		xPoints[i] = - dG(xPoints[i]) / ddG(xPoints[i])


remezError  = np.polyval(P,x)
plt.plot(x, remezError  - f(x), 'b-')
for i in range(len(xPoints)) :
	plt.plot([xPoints[i], xPoints[i]], [-eps, eps], 'y-*')


L = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
P = P[::-1]
for i in range(len(P)) :
	print("#define COEFF_P"+str(N)+"_"+str(L[i])+" "+str(P[i]))




	
plt.show()