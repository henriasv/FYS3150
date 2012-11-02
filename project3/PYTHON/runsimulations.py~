import subprocess
# input format for program: $ <exec> <N> <range> <method> <output> <N_cores>


exec_path = "/mn/felt/u8/henriasv/Dropbox/Henrik/Emner/FYS3150/FYS3150-code/project3/dist/Debug/GNU-Linux-x86/project3"
outpath = "/home/scratch/henriasv/FYS3150/project3data"
N_threads = "8"

# direct methods:

methods_dir = ["gaulag", "gauleg"]
N_dir = []# ["3","5","8"]#["%d" %i for i in range(10, 35, 5)]
print N_dir;
ranges_dir = ["3"]

for method in methods_dir:
	for N in N_dir:
		for ranges in ranges_dir:
			subprocess.call([exec_path, N, ranges, method, outpath, N_threads])
			print exec_path, N, ranges, method, outpath, N_threads


# monte carlo methods:

methods_mc = ["bf_mc", "is_mc"]
N_mc = ["%d" % 10**i for i in range(2, 3)];
print N_mc
ranges_mc = ["3"]

for method in methods_mc:
	for N in N_mc:
		for ranges in ranges_mc:
			subprocess.call([exec_path, N, ranges, method, outpath, N_threads])
			print exec_path, N, ranges, method, outpath, N_threads

