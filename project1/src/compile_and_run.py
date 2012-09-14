import subprocess
import sys
import os


filename = os.getcwd()+"/"+sys.argv[1];
codeFile = filename + ".cpp";
libFile = os.getcwd()+"/"+"lib.cpp";

a = subprocess.call(['c++', '-c', "-Wall", codeFile, libFile]);
if not a==0:
	print "Could not compile, quitting"
	sys.exit(1)
subprocess.call(['c++', '-o', filename, filename+'.o', "lib.o"])
subprocess.call([filename])
