#!/usr/bin/python
import os
import string
import sys
from time import *

def IntToBool(Integer):
	if Integer==0 :
		return "false"
	else :
		return "true"

###################################################################################################
def set_param(dev, shared_mem, doub_prec, mat_low, mat_step, mat_up, average_runs, lu ):

	prm = open("./prm/step-1-Decomp.prm","r")
	prm_list = prm.readlines()
	prm.close()
	i = 0
	for line in prm_list:	
		if(str.find(line,"set Device") != -1):
			prm_list[i] = "  set Device        = "+str(dev)+"\n"
		if(str.find(line,"set Shared Memory") != -1):
			prm_list[i] = "  set Shared Memory = "+IntToBool(shared_mem)+"\n"
		if(str.find(line,"set Double-Precision") != -1):
			prm_list[i] = "  set Double-Precision          = "+IntToBool(doub_prec)+"\n"
		if(str.find(line,"set Matrix size - lower limit") != -1):
			prm_list[i] = "  set Matrix size - lower limit = "+str(mat_low)+"\n"
		if(str.find(line,"set Matrix size - step size") != -1):
			prm_list[i] = "  set Matrix size - step size   = "+str(mat_step)+"\n"
		if(str.find(line,"set Matrix size - upper limit") != -1):
			prm_list[i] = "  set Matrix size - upper limit = "+str(mat_up)+"\n"
		if(str.find(line,"set LU decomposition") != -1):
			prm_list[i] = "  set LU decomposition          = "+IntToBool(lu)+"\n"

		i +=1

	prm = open("./prm/step-1-Decomp.prm","w")
	prm.writelines(prm_list)
	prm.close()




#DEVICE	= 0
#SHARED_MEMORY=0
#DOUBLE_PRECISION=0
mat_low=256
mat_step=64
mat_up=513
dev=0
shmem=0
doub_prec=1
lu=1
average_runs=3

eingabe="K"

while (eingabe.upper()!="E"):
	print("-------------------------------------------------------")
	print("Choose which parameter to change")	
	print("-------------------------------------------------------\n\n")
	print("(D)evice\t\t\t\t set to: %i" % dev)
	print("(S)hared memory\t\t\t\t set to: %i" % shmem)
	print("Double (p)recision\t\t\t set to: %i" % doub_prec)
	print("(L)U Decomposition\t\t\t set to: %i" % lu)
	print("Matrix l(o)wer limit\t\t\t set to: %i" % mat_low)
	print("Matrix (u)pper limit\t\t\t set to: %i" % mat_up)
	print("Matrix s(t)ep size\t\t\t set to: %i" % mat_step)
	print("(A)veraging runs\t\t\t set to: %i" % average_runs)
	print("T(r)oubleshooting")
	print("Save and (E)xit")

	print("\n-------------------------------------------------------\n")
	eingabe = raw_input("Choice: ")
	
	if (eingabe.upper() == "D"):
		dev = int(raw_input("New value: "))
	elif (eingabe.upper() == "S"):
		shmem = int(raw_input("New value: "))
	elif (eingabe.upper() == "P"):
		doub_prec = int(raw_input("New value: "))
	elif (eingabe.upper() == "L"):
		lu = int(raw_input("New value: "))
	elif (eingabe.upper() == "O"):
		mat_low = int(raw_input("New value: "))
	elif (eingabe.upper() == "U"):
		mat_up = int(raw_input("New value: "))
	elif (eingabe.upper() == "T"):
		mat_step = int(raw_input("New value: "))
	elif (eingabe.upper() == "A"):
		average_runs = int(raw_input("New value: "))
	elif (eingabe.upper() == "R"):
		print("1.The script doesn't change my parameter file!\n>>>The script needs the parameter file to be located in a folder named 'prm'. Additionally the file itself has to be named 'step-1-Decomp.prm' and it must exist!")
		print("2. My program's parameters differ from the ones I set in the script!\n>>>If not further prepared the script has to be located in the launch directory of the program. This might be the build-desktop-directory")
		print("3. The script and/or the program crashes\n>>>Be sure to enter only valid values for the corresponding variables. This might be 0 and 1 (for bools) or general integers.")
	print("\n")

set_param(dev,shmem,doub_prec,mat_low,mat_step,mat_up,average_runs,lu)









