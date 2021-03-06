from scitools.easyviz.matplotlib_ import plot, hold, figure, loglog, xlabel, ylabel, title, legend
import numpy as np
import os

def analytic_solution(x):
	return 1-x+x*np.exp(-10)-np.exp(-10*x)

def plotSolution(path, error):
	"""
	Definition space of the solution is assumed to be [0, 1)
	Length of loaded array is assumed to be N = N_sim+2, where 2 comes 
	from adding the boundaries after calculations.
	h = 1/(N_sim+1) = 1/(N-1)
	Plotting only solutions with high steplength to avoid crash of gnuplot
	"""
	dataSolution = np.fromfile(path, dtype=np.float64);
	N = len(dataSolution);
	h = 1.0/(N-1);
	x = np.linspace(0, 1, N);
	solution = analytic_solution(x);
	if N<10000:
		figure()
		plot(x, dataSolution, 'o');
		hold('on');
		plot(x, solution);
		title('Solution to -u\'\'=f(x), h=%6.2g' %h)
		xlabel('x')
		ylabel('u(x)')
		legend('Approximation', 'Analytic solution')
	
	error_array = solution-dataSolution
	error_array = np.abs(error_array/solution)
	error_array[0] = 0;
	error_array[len(error_array)-1] =0;
	
	error.append([h, np.max(error_array)])

path = "/scratch/henriasv/FYS3150/project1/"
error = [];
for name in os.listdir(path):
	if not name.endswith('.txt'):
		plotSolution(path+name, error)
error = np.asarray(error)
error = np.log10(error)
figure()
plot(error[:, 0], error[:, 1], '*')
title('maximum relative error with different step lengths')
xlabel('log(h)')
ylabel('log($\epsilon$)')


raw_input("press enter")
