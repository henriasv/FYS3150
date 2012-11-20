# Script running ising code in parallell. Each N value serial, but T values in each N are in parallell

import subprocess
from threading import Thread, Lock
from numpy import linspace
from scipy.stats import norm
import sys

if not (raw_input("sure?") == "yes"):
	sys.exit(0)

path = "/mn/felt/u8/henriasv/Dropbox/Henrik/Emner/FYS3150/FYS3150-code/IsingModel2D/dist/Debug/GNU-Linux-x86/isingmodel2d";

n_mcs = 5*100000
N_values = [2, 20, 30, 40, 60, 80, 100];
T_values = 0.12*norm.ppf(linspace(0.01, 0.99, 20)) + 2.269

def simulationsT(N_values, T):
	outputs = [];
	for N in N_values:
		print "N", N, "T", T
		proc = subprocess.Popen([path, "%d"%N , "%f" % T, "%d" % n_mcs], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
		output, unused_err = proc.communicate()  # buffers the output
		retcode = proc.poll()                    # ensures subprocess termination
		outputs.append(output)
	for output in outputs:
		print output
processes = []

for T in T_values:
	processes.append(Thread(target=simulationsT, args=(N_values, T)))
	processes[-1].start()
for p in processes: p.join()
print('done.')
