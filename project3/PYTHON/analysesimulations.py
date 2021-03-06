import os
import numpy as np
from scitools.easyviz.matplotlib_ import plot, loglog, legend, hold, title, xlabel, ylabel, figure, semilogx
	

class Simulation:
	def __init__(self, path, folder):
		self.exact_solution = 5*np.pi**2/(16**2);
		self.path = path
		self.folder = folder
		self.N = 0
		self.N_threads = 0
		self.energy = 0
		self.error = 0
		self.resolution =  0
		self.method = " "
		self.getparameters()


	def getparameters(self):
		strings = self.folder.split("#")
		for substring in strings:
			tmpStr = substring.split("%")
			try:
				if tmpStr[0] == "N":
					self.N = float(tmpStr[1]);
					self.resolution = self.N;
			except IndexError:
				print tmpStr;
		for substring in strings:
			tmpStr = substring.split("%")
			try:
				if tmpStr[0] == "method":
					if tmpStr[1] == "is_mc":
						self.ID = "is_mc"
						self.method = "Importance sampling Monte Carlo"
						self.resolution = self.N**(1.0/6)
					if tmpStr[1] == "bf_mc":
						self.ID = "bf_mc"
						self.method = "Brute force Monte Carlo"
						self.resolution = self.N**(1.0/6)
					if tmpStr[1] == "gauleg":
						self.ID = "gauleg"
						self.method = "Gauss-Legendre Quadrature"
						self.N = self.N**6
					if tmpStr[1] == "gaulag":
						self.ID = "gaulag" 
						self.method = "Gauss-Laguerre Quadrature"
						self.N = self.N**6
				if tmpStr[0] == "N_threads":
					self.N_threads = float(tmpStr[1]);
			except IndexError:
				print tmpStr;
		data = np.fromfile(self.path + "/" + self.folder + "/output.bin")
		self.energy = data[0]
		self.cpu_time = data[1]
		self.error = abs(self.energy-self.exact_solution)/exact_solution



	def __str__(self):
			return "-"*40 + "\n" +  "Simulation with " + self.method + "\n" + \
				"E = " + "%g\n" %self.energy + "Error = %g\n" %self.error + \
				"cpu-time = " + "%g" % self.cpu_time + "\n" + \
				"Threads = " + "%d" % self.N_threads + "\n" +\
				"N = " + "%g"%self.N + "\n" + "Resolution = " + "%g" %self.resolution  \
				 + "\n" + "-"*40

def sort_data(data):
	ordered = False
	while not ordered:
		for i in range(len(data)-1):
			if data[i, 0] <= data[i+1, 0]:
				ordered = True
			else:
				tmpdata0 = data[i, 0]
				tmpdata1 = data[i, 1]
				data[i, 0] = data[i+1, 0]
				data[i, 1] = data[i+1, 1]
				data[i+1, 0] = tmpdata0
				data[i+1, 1] = tmpdata1
				ordered = False
				break
		

if __name__ == '__main__':

	path = "/home/scratch/henriasv/FYS3150/project3data"
	paths = os.listdir(path)
	exact_solution = 5*np.pi**2/(16**2);
	print "Exact solution %g" %exact_solution

	simulations = [];
	for folder in paths:		
		tmp = Simulation(path, folder)
		print tmp
		simulations.append(tmp)


	################################################################
	# Plots: 	Precision as a function of number of function evaluations for all methods
	# 			cpu-time as a function of function evaluations for all methods
	# 			Best energies (table)
	#

	
	# Plotting function evaluations and precision
	prec_gaulag = []
	prec_gauleg = []
	prec_is_mc  = []
	prec_bf_mc  = []

	for sim in simulations:
		if sim.ID == "gaulag":
			prec_gaulag.append([sim.N, sim.error])
		if sim.ID == "gauleg":
			prec_gauleg.append([sim.N, sim.error])
		if sim.ID == "bf_mc":
			prec_bf_mc.append([sim.N, sim.error])
		if sim.ID == "is_mc":
			prec_is_mc.append([sim.N, sim.error])
	prec_gaulag = np.asarray(prec_gaulag)
	prec_gauleg = np.asarray(prec_gauleg)
	prec_is_mc = np.asarray(prec_is_mc)
	prec_bf_mc = np.asarray(prec_bf_mc)
	sort_data(prec_gaulag)
	sort_data(prec_gauleg)
	sort_data(prec_is_mc)
	sort_data(prec_bf_mc)

	figure()
	loglog(prec_gaulag[:,0], prec_gaulag[:,1])
	hold("all")
	loglog(prec_gauleg[:,0], prec_gauleg[:,1])
	loglog(prec_bf_mc[:,0], prec_bf_mc[:,1])
	loglog(prec_is_mc[:,0], prec_is_mc[:,1])

	legend("Gauss Laguerre", "Gauss Legendre", "Brute force Monte Carlo", "Importance sampling Monte Carlo")
	xlabel("Number of function evaluations")
	ylabel("Relative error")
	title("Relative error vs number of calculations for different methods")
	hold("off")

	# Plotting function evaluations and cpu_time
	prec_gaulag = []
	prec_gauleg = []
	prec_is_mc  = []
	prec_bf_mc  = []

	for sim in simulations:
		if sim.ID == "gaulag":
			prec_gaulag.append([sim.N, sim.cpu_time])
		if sim.ID == "gauleg":
			prec_gauleg.append([sim.N, sim.cpu_time])
		if sim.ID == "bf_mc":
			prec_bf_mc.append([sim.N, sim.cpu_time])
		if sim.ID == "is_mc":
			prec_is_mc.append([sim.N, sim.cpu_time])
	prec_gaulag = np.asarray(prec_gaulag)
	prec_gauleg = np.asarray(prec_gauleg)
	prec_is_mc = np.asarray(prec_is_mc)
	prec_bf_mc = np.asarray(prec_bf_mc)
	sort_data(prec_gaulag)
	sort_data(prec_gauleg)
	sort_data(prec_is_mc)
	sort_data(prec_bf_mc)

	figure()
	loglog(prec_gaulag[:,0], prec_gaulag[:,1])
	hold("all")
	loglog(prec_gauleg[:,0], prec_gauleg[:,1])
	loglog(prec_bf_mc[:,0], prec_bf_mc[:,1])
	loglog(prec_is_mc[:,0], prec_is_mc[:,1])
	legend("Gauss Laguerre", "Gauss Legendre", "Brute force Monte Carlo", "Importance sampling Monte Carlo")
	xlabel("Number of function evaluations")
	ylabel("cpu-time")
	hold("off")
	

	# Plotting cpu_time and precision
	prec_gaulag = []
	prec_gauleg = []
	prec_is_mc  = []
	prec_bf_mc  = []

	for sim in simulations:
		if sim.ID == "gaulag":
			prec_gaulag.append([sim.cpu_time, sim.error])
		if sim.ID == "gauleg":
			prec_gauleg.append([sim.cpu_time, sim.error])
		if sim.ID == "bf_mc":
			prec_bf_mc.append([sim.cpu_time, sim.error])
		if sim.ID == "is_mc":
			prec_is_mc.append([sim.cpu_time, sim.error])
	prec_gaulag = np.asarray(prec_gaulag)
	prec_gauleg = np.asarray(prec_gauleg)
	prec_is_mc = np.asarray(prec_is_mc)
	prec_bf_mc = np.asarray(prec_bf_mc)
	sort_data(prec_gaulag)
	sort_data(prec_gauleg)
	sort_data(prec_is_mc)
	sort_data(prec_bf_mc)

	figure()
	loglog(prec_gaulag[:,0], prec_gaulag[:,1])
	hold("all")
	loglog(prec_gauleg[:,0], prec_gauleg[:,1])
	loglog(prec_bf_mc[:,0], prec_bf_mc[:,1])
	loglog(prec_is_mc[:,0], prec_is_mc[:,1])
	legend("Gauss Laguerre", "Gauss Legendre", "Brute force Monte Carlo", "Importance sampling Monte Carlo")
	xlabel("cpu-time")
	ylabel("Relative error")
	hold("off");


	# Actual convergence
	

	prec_gaulag = []
	prec_gauleg = []
	prec_is_mc  = []
	prec_bf_mc  = []

	for sim in simulations:
		if sim.ID == "gaulag":
			prec_gaulag.append([sim.N, sim.energy])
		if sim.ID == "gauleg":
			prec_gauleg.append([sim.N, sim.energy])
		if sim.ID == "bf_mc":
			prec_bf_mc.append([sim.N, sim.energy])
		if sim.ID == "is_mc":
			prec_is_mc.append([sim.N, sim.energy])
		exact = sim.exact_solution;
	prec_gaulag = np.asarray(prec_gaulag)
	prec_gauleg = np.asarray(prec_gauleg)
	prec_is_mc = np.asarray(prec_is_mc)
	prec_bf_mc = np.asarray(prec_bf_mc)
	sort_data(prec_gaulag)
	sort_data(prec_gauleg)
	sort_data(prec_is_mc)
	sort_data(prec_bf_mc)

	figure()
	semilogx(prec_gaulag[:,0], prec_gaulag[:,1])
	hold("all")
	semilogx(prec_gauleg[:,0], prec_gauleg[:,1])
	semilogx(prec_bf_mc[:,0], prec_bf_mc[:,1])
	semilogx(prec_is_mc[:,0], prec_is_mc[:,1])
	semilogx([1e2, 1e10], [exact, exact], '--')
	legend("Gauss Laguerre", "Gauss Legendre", "Brute force Monte Carlo", "Importance sampling Monte Carlo", "exact")
	xlabel("Function evaluations")
	ylabel("Calculated energy")
	title("Convergence of the different methods")
	hold("off");
	
	print "Results with highest resolution for all methods"
	print "Gauss Legendre: %.8f" % prec_gauleg[len(prec_gauleg)-1, 1], "Evaluations: %g" % prec_gauleg[len(prec_gauleg)-1, 0]
	print "Gauss Laguerre: %.8f" % prec_gaulag[len(prec_gaulag)-1, 1], "Evaluations: %g" % prec_gaulag[len(prec_gaulag)-1, 0]
	print "Brute force Monte Carlo: %.8f" % prec_bf_mc[len(prec_bf_mc)-1, 1], "Evaluations: %g" % prec_bf_mc[len(prec_bf_mc)-1, 0]
	print "Importance sampling Monte Carlo %.8f" % prec_is_mc[len(prec_is_mc)-1, 1], "Evaluations: %g" % prec_is_mc[len(prec_is_mc)-1, 0]
	print "Exact: ", "%.8f" % exact

	raw_input("press enter")



































