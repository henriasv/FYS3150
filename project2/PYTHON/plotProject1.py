from scitools.easyviz.matplotlib_ import plot, hold, figure, loglog, xlabel, ylabel, title, legend, semilogy, semilogx, closefigs
from scipy.integrate import trapz
from matplotlib import rc
import os
import numpy as np
import sys

path = "/scratch/henriasv/FYS3150/project2/"  + sys.argv[1] + "/";

class Simulation:
	def __init__(self, path, max_states=3):
		self.eigenvalues=None;
		self.rho = None;
		self.states = [[]]*max_states;
		self.omega = None;
		self.N = None;
		self.names = [[]]*max_states;
		print self.states
		

		for name in os.listdir(path):
			if "grid" in name:
				self.rho = np.fromfile(path+name, dtype=np.float64)
			if "eigenvalues" in name:
				self.eigenvalues = np.fromfile(path+name, dtype=np.float64)
			if "results" in name:
				infile = open(path+name);
				for line in infile:
					content = line.split("=");
					if "omega" in content[0].strip():
						self.omega = float(content[1]);

		for name in os.listdir(path):
			if "state" in name:
				state_number = int(name.strip(".bin").strip("state"))
				state = np.fromfile(path+name, dtype=np.float64)/self.rho;
				
				# Normalizing wave function
				state = state/np.sqrt(trapz(state**2, self.rho));
				self.states[state_number] = state;
				
				# setting number of grid points for reference
				self.N = len(state)
		plot(self.rho, self.states[1]);


class StatePlotter:
	def __init__(self, simulations):
		self.simulations = simulations
		min_states = False;
		for sim in simulations:
			if not min_states:
				min_states = len(sim.states);
			elif min_states>len(sim.states):
				min_states = len(sim.states);
			else:
				None

	def plotByOmega(self, omegalimit = 0.1):
		for sim in self.simulations:
			for i in range(len(sim.states)):
				if sim.omega >= omegalimit:
					figure(i);
					plot(sim.rho, sim.states[i]**2, legend="$\omega= %g, \lambda= %g$" % (sim.omega, sim.eigenvalues[i]), xlabel=r'$\frac{1}{\alpha}r$', \
							ylabel=r'$\alpha|\psi_n|^2$')
					title(r"$\psi_{%d}$" % i);
					hold('on')
		raw_input("press enter")
		closefigs();

	def plotBySimulation(self):
		for i in range(len(self.simulations)):
			sim = self.simulations[i];
			figure(i);
			for j in range(len(sim.states)):
				plot(sim.rho, sim.states[j]**2, legend=r"$\psi_{%d}, \lambda= %g$" % (j, sim.eigenvalues[j]), xlabel=r'$\frac{1}{\alpha}r$', ylabel=r'$\alpha|\psi_n|^2$')
				title("$\omega = %g$" %sim.omega)
				hold('on')
		raw_input("press enter")
		closefigs()

class Simulations:
	def __init__(self, path):
		self.sims = []
		for name in os.listdir(path):
			self.sims.append(Simulation(path+name+"/"))	

if __name__ == "__main__":
	sims = None
	if sys.argv[2] == "multiple":
		tmp = Simulations(path);
		sims = tmp.sims
	if sys.argv[2] == "single":
		sims = [Simulation(path)]
	else:
		print "Supply mode: single/multiple";
	
	plotter = StatePlotter(sims);
	plotter.plotByOmega();
	plotter.plotBySimulation();
