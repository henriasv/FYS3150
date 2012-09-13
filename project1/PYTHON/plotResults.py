from scitools.easyviz.matplotlib_ import plot, hold, figure
import numpy as np
import os

def analytic_solution(x):
	return 1-x+x*np.exp(-10)-np.exp(-10*x)

def plotSolution(path):
	dataSolution = np.fromfile(path, dtype=np.float64);
	N = len(dataSolution);
	h = 1.0/(N+1);
	x = np.linspace(h, 1-h, N);
	solution = analytic_solution(x);
	figure();
	plot(x, dataSolution, 'o');
	hold('on');
	plot(x, solution);

path = "/scratch/henriasv/FYS3150/project1/"
for name in os.listdir(path):
	plotSolution(path+name)


raw_input("press enter")
