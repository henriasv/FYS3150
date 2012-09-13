import subprocess
import sys
import os


filename = os.getcwd()+"/"+sys.argv[1];
codeFile = filename + ".cpp";

a = subprocess.call(['c++', '-c', codeFile]);
if not a==0:
	print "Could not compile, quitting"
	sys.exit(1)
subprocess.call(['c++', '-o', filename, filename+'.o'])
subprocess.call([filename])
