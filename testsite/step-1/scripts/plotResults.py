#!/usr/bin/python
import os
import string
import sys
from time import *


def plot_results(output, xlabel, ylabel, xcol, ycol, files, ytitle):
	plot = "set terminal " + TERMINAL \
	+ "\nset output \'" + output + "\'" \
	+ "\nset xlabel \'" + xlabel + "\'" \
	+ "\nset ylabel \'" + ylabel + "\'" \
	+ "\n" + GLOBAL_FORMAT  \
	+ "\nplot " 
	for i in range(0,len(files)-1):
		plot = plot + "\'" + files[i] + "\' u " + str(xcol) + ":" + str(ycol) + " " + FORMAT + " title \'" + ytitle[i] + "\', "

	plot = plot + "\'" + files[len(files)-1] + "\' u " + str(xcol) + ":" + str(ycol) + " " + FORMAT + " title \'" + ytitle[len(files)-1] + "\'"
	
	os.system("echo \"" + plot + "\"")
	os.system("echo \"" + plot + "\" | gnuplot -persist")


if len(sys.argv)==1:
	print "Choose one argument to be plotted: device, precision, memory"
else:
	PLOTFAC = sys.argv[1]
	
	RPATH = "test_me"

	TERMINAL = "wxt"

	GLOBAL_FORMAT = ""
	FORMAT = "with linespoints lw 1"
	XCOL = 1
	YCOL = 3

	dev="dev0"
	mem="shmem"
	prec="float"

	FILENAME = []
	FILENAME.append("./"+ RPATH + "/chol_fac_times_" + dev + "_" + mem + "_" + prec + ".dat")

	if PLOTFAC == "device":
		dev="dev1"
		PLOTCAP = ["Device 0", "Device 1"]
	if PLOTFAC=="memory":
		mem="noshmem"
		PLOTCAP = ["Shared Memory", "No Shared Memory"]
	if PLOTFAC=="precision":
		prec="double"
		PLOTCAP = ["Float Precision" , "Double Precision"]

	FILENAME.append("./"+ RPATH + "/chol_fac_times_" + dev + "_" + mem + "_" + prec + ".dat")

	plot_results("./"+ RPATH + "/out1", "matrix size", "computation time", XCOL, YCOL, FILENAME, PLOTCAP)


#os.system("echo \"" + plot + "\"")
#os.system("echo \"" + plot + "\" | gnuplot -persist"  )
	

