#!/usr/bin/python
import os
import string
from time import *

TERMINAL = "wxt"
FORMAT = "with linespoints lw 1"
FILENAME = "mergedResults.dat"
XCOL = 1
YCOL = 2

def plot_results(output, xlabel, ylabel, xcol, ycol, ytitle,files):
	plot = "set terminal " + TERMINAL \
	+ "\nset output \'" + output + "\'" \
	+ "\nset xlabel \'" + xlabel + "\'" \
	+ "\nset ylabel \'" + ylabel + "\'" \
	+ "\nplot " 
	for i in range(0,len(ycol)-1):
		plot = plot + "\'" + FILENAME + "\' u " + str(XCOL) + ":" + str(YCOL) + " " + FORMAT + " title \'" + ytitle[i] + "\', "

	plot = plot + "\'" + FILENAME + "\' u " + str(xcol) + ":" + str(ycol[len(ycol)-1]) + " " + FORMAT + " title \'" + ytitle[len(ycol)-1] + "\'"
	
	os.system("echo \"" + plot + "\"")
	os.system("echo \"" + plot + "\" | gnuplot -persist")

plot_results("out1", "matrix size", "computation time", 1, [2, 3], ["erste Spalte", "zweite Spalte"], ["chol_fac_times_0_0_0.dat", "chol_fac_times_0_0_1.dat"] )
#os.system("echo \"" + plot + "\"")
#os.system("echo \"" + plot + "\" | gnuplot -persist"  )
	

