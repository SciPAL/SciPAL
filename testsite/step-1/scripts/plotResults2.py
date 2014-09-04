#!/usr/bin/python
import os
import string
import sys
from time import *


prm = open("./prm/step-1-Decomp.prm","r")
prm_list = prm.readlines()
prm.close()
for line in prm_list:	
		if(str.find(line,"set Run directory") != -1):
			PATH=(line.split(" "))[len(line.split(" "))-1]
			PATH=PATH.split("\n")[0]

output_files=[]


def KonfigurationsPlain(dev,mem,prec,lu):
	if mem==0:
		memstr="no shared memory"
	else:
		memstr="shared memory"
	if prec==0:
		precstr="float-precision"
	else:
		precstr="double-precision"
	if dev==0:
		devstr="Device 0"
	else: 
		devstr="Device 1"
	if lu==0:
		methstr="Cholesky decomposition"
	else:
		methstr="LU decomposition"

	return memstr + ",\t" + precstr + ",\t" + devstr + ",\t" + methstr 



def Filename(dev, mem, prec, lu):
	if mem==0:
		memstr="noshmem"
	else:
		memstr="shmem"
	if prec==0:
		precstr="float"
	else:
		precstr="double"
	if lu==0:
		methstr="chol"
	else:
		methstr="lu"

	out=PATH +"/"+ methstr + "_fac_times_dev" + str(dev) + "_" + memstr + "_" + precstr + ".dat_avg"

	return out

def plot_results(output, xlabel, ylabel, xcol, ycol, files, ytitle):
	plot = "set terminal " + TERMINAL \
	+ "\nset output \'" + output + "\'" \
	+ "\nset xlabel \'" + xlabel + "\'" \
	+ "\nset ylabel \'" + ylabel + "\'" \
	+ "\n" + GLOBAL_FORMAT  \
	+ "\nplot " 
	for i in range(0,len(files)-1):
		plot = plot + "\'" + files[i] + "\' u " + str(xcol) + ":" + str(ycol) + FORMAT + " title \'" + ytitle[i] + "\', "
		plot = plot + "\'" + files[i] + "\' u " + str(xcol) + ":" + str(ycol) + ":" + str(ycol+1) + " w yerrorbars notitle , " 
		

	plot = plot + "\'" + files[len(files)-1] + "\' u " + str(xcol) + ":" + str(ycol) + " " + FORMAT + " title \'" + ytitle[len(files)-1] + "\', "
	plot = plot + "\'" + files[len(files)-1] + "\' u " + str(xcol) + ":" + str(ycol) + ":" + str(ycol+1) + " w yerrorbars notitle " 
	plot = plot + "\nset output"
	os.system("echo \"" + plot + "\"")
	os.system("echo \"" + plot + "\" | gnuplot -persist")


if len(sys.argv)==1:
	print "Choose one argument to be plotted: device, precision, memory, decompType"
else:
	PLOTFAC = sys.argv[1]
	
	RPATH = "test_me"

	TERMINAL = "postscript"

	GLOBAL_FORMAT = ""
	FORMAT = " with linespoints lw 1"
	XCOL = 1
	
	YCOL = [3,5,7,9,11,13,15]
	YLABEL = ["Cpu factorization time", "Pure GPU factorization time", "GPU fact. time incl data transfer", "diag\_update", "factorize\_diagonal\_block", "lo\_update", "strip\_update"]

	for i in range(0, len(YCOL)):

		output = PLOTFAC + str(YCOL[i]) + ".ps"
		
		dev=0
		mem=0
		prec=0
		lu=0
		configs = []
		configs.append(KonfigurationsPlain(dev,mem,prec,lu));
		FILENAME = []
		FILENAME.append(Filename(dev,mem,prec,lu))
	
		if PLOTFAC == "device":
			dev=1
			PLOTCAP = ["Device 0", "Device 1"]
		if PLOTFAC=="memory":
			mem=1
			PLOTCAP = ["Shared Memory", "No Shared Memory"]
		if PLOTFAC=="precision":
			prec=1
			PLOTCAP = ["Float Precision" , "Double Precision"]
		if PLOTFAC == "decompType":
			lu=1
			PLOTCAP = ["Cholesky decomposition","LU Decomposition"]
			
	
		FILENAME.append(Filename(dev,mem,prec,lu))
		configs.append(KonfigurationsPlain(dev,mem,prec,lu));
		plot_results(PATH + "/"+output, "Matrix Size", YLABEL[i], XCOL, YCOL[i], FILENAME, PLOTCAP)
		output_files.append(output)

#os.system("echo \"" + plot + "\"")
#os.system("echo \"" + plot + "\" | gnuplot -persist"  )
	
text_lst = []
text_lst.append("\documentclass[10pt]{article}\n")
text_lst.append("\usepackage{graphicx}\n")
text_lst.append("\usepackage{float}\n")
text_lst.append("\usepackage{subfig}\n")
text_lst.append("\graphicspath{{"+PATH+"/}}\n")
text_lst.append("\\begin{document}\n")
text_lst.append("\centering\\textbf{Results of step-1 -- Parameter: "+str(PLOTFAC)+"}\n")
subfig=1
for outfile in output_files:
	if (subfig % 2):
		text_lst.append("\\begin{figure}[H]\n\t\\begin{center}\n")
	text_lst.append("\t\t\subfloat{\includegraphics[angle=-90,width=0.45\\textwidth]{"+RPATH+"/"+outfile+"}}\n")
	if (subfig % 2 == 0):
		text_lst.append("\t\end{center}\n\end{figure}\n")
	subfig+=1
if (subfig % 2==0):
	text_lst.append("\t\end{center}\n\end{figure}\n")
#Tabelle mit Ergebnissen
#text_lst.append("\clearpage")
for filenr in range(0,len(FILENAME)):
	data = open(FILENAME[filenr],"r")
	data_lines = data.readlines()
	data.close()	
	text_lst.append("\\textsc{"+configs[filenr]+"}\n")
	text_lst.append("\\begin{table}[H]\n\t\centering\n\t\\begin{tabular}{|c|c|c|c|}\n")
	text_lst.append("\t\t\hline\n \t\tMatrix size & CPU fac time [s] & Pure GPU fac time [s] & GPU fac t incl data trans. [s]\\\\\n\t\t\hline\n\t\t")
	for i in data_lines:
		line = i.split("\t")
		for value in range(0,len(line)/2-4):
			if (value):
				text_lst.append(line[value*2]+" $\pm$ ")
				text_lst.append(line[value*2+1])
			else:
				text_lst.append(line[value*2])
			if (value < len(line)/2-5):
				text_lst[len(text_lst)-1]+= " & "
		text_lst.append("\\\\\hline\n\t\t")
	text_lst.append("\\end{tabular}\n\\end{table}")
	text_lst.append("\\begin{table}[H]\n\t\centering\n\t\\begin{tabular}{|c|c|c|c|}\n")
	text_lst.append("\t\t\hline\n \t\tdiag\_update [ms] & factorize\_diag\_block [ms] & lo\_update [ms] & strip\_update [ms]\\\\\n\t\t\hline\n\t\t")
	for i in data_lines:
		line = i.split("\t")
		for value in range(4,len(line)/2):
			text_lst.append(line[value*2]+" $\pm$ ")
			text_lst.append(line[value*2+1])
			if (value < len(line)/2-1):
				text_lst[len(text_lst)-1]+= " & "
		text_lst.append("\\\\\hline\n\t\t")
	text_lst.append("\\end{tabular}\n\\end{table}")
text_lst.append("\end{document}\n")



latex = open("results.tex","w")
latex.writelines(text_lst)
latex.close()
os.system("latex results.tex")






