#!/usr/bin/env python
 
#import random_pkg
import numpy as np
from scipy.stats import poisson

# run parameters
#mutn_rate = 0.01
mutn_rate = 1.2e-6
ngen = 100
pop_size = 1000

# seed the random number generators
#rns = int(raw_input("enter the random number seed:"))
for rns in range(1,10):
	np.random.seed(rns)
	#random_pkg.random_pkg.randomnum(-rns)

	cumR = 0
	cumS = 0
	S = []


	for i in range(ngen):
		#random_pkg.random_pkg.randomnum(-rns-i)
		#for j in range(pop_size): 
		#	S += [ random_pkg.random_pkg.poisson(mutn_rate) ]
		R = poisson.rvs(mutn_rate, size=pop_size)
		cumR += sum(R)
		#cumS += sum(S)
		#del S[:]

	#print 'RNS:', rns,'Mersenne Twister:', cumR, 'Numerical Recipes:', cumS
	print 'RNS:', rns,'Mersenne Twister:', cumR

