from numpy import asarray, sinh, cosh, exp

def expE(T):
	return -8*sinh(8./T)/(cosh(8./T)+3)/4

def expM(T):
	return (2*exp(8./T) + 4)/(cosh(8./T)+3)/4

T = asarray([1, 2.4])

E = expE(T)
M = expM(T)

for i in range(len(T)):
	print "--------------------------"
	print "Expectation values per spin for L = 2 and T = %g" % T[i]
	print "<E> = %g" % E[i]
	print "<|M|> = %g" % M[i]
