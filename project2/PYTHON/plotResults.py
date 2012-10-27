from scitools.easyviz.matplotlib_ import plot, hold, figure, loglog, xlabel, ylabel, title, legend, semilogy
from matplotlib import rc
from scipy.integrate import trapz
import numpy as np
import os
import sys

rc('text', usetex=True)
rc('font', family='serif')


def plotSolution(path, name, rho, eigenvalues):
	dataSolution = np.fromfile(path+name, dtype=np.float64);
	N = len(dataSolution);
	x = rho;
	dataSolution = dataSolution/np.sqrt(trapz(dataSolution**2/x**2, x))
	state_number = round(float((name.strip(".bin").strip("state"))));
	integral = trapz(dataSolution**2/x**2, x)
	print state_number;
	if N<10000:
		figure()
		plot(x, dataSolution**2/x**2);
		hold('on');
		xlabel(r'$\rho$')
		ylabel(r'$R^2(\rho)$')	
		title(name.strip('.bin') + ", " + '$\lambda =$ %g' % eigenvalues[state_number])
		legend("integral = %g" % integral)
	else:
		print "Stopped plot, N>10000";

path = "/scratch/henriasv/FYS3150/project2/"  + sys.argv[1];

eigenvalues = None
rho = None;
sqsum = None;
for name in os.listdir(path):
	if "eigenvalues" in name:
		eigenvalues = np.fromfile(path+name, dtype=np.float64)
		print eigenvalues

	if "grid" in name:
		print path+name
		rho = np.fromfile(path+name, dtype=np.float64);
		print rho

	if "sqsum" in name:
		sqsum = [];
		fileInput = open(path + name)
		for line in fileInput:
			sqsum.append(float(eval(line)))
		sqsum = np.asarray(sqsum)
		figure()
		semilogy(sqsum)
		title("Sum of squared off-diagonal elements")
	

for name in os.listdir(path):
	if "state" in name:
		plotSolution(path, name, rho, eigenvalues)
	


raw_input("press enter")
