#!/usr/bin/python
import os
import string
import sys
from time import *




#DEVICE	= 0
#SHARED_MEMORY=0
#DOUBLE_PRECISION=0
MATRIX_LOW=256
MATRIX_STEP=64
MATRIX_UP=513
AVERAGE_RUNS=int(sys.argv[1])

PATH="./test_me/"
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

def IntToBool(Integer):
	if Integer==0 :
		return "false"
	else :
		return "true"

def Filename(dev, mem, prec,lu):
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

	out=PATH + methstr+"_fac_times_dev" + str(dev) + "_" + memstr + "_" + precstr + ".dat"

	return out

for DEVICE in range(0,2):
	for SHARED_MEMORY in range(0,2):
		for DOUBLE_PRECISION in range(0,2):
			for LU in range(0,2):
			
#				act_run = [AVERAGE_RUNS],[],[]]
				arr=[]
				for run in range(0,AVERAGE_RUNS):
					print "---------------------------------------------------------------------------------------\n"
					print "\t\t\t\tKonfiguration: "+str(LU+DOUBLE_PRECISION*2+SHARED_MEMORY*4+DEVICE*8+1) + " von 16" 
					print "\t\t\t\tMitteln: "+str(run+1)+ " von " + str(AVERAGE_RUNS)
					print "---------------------------------------------------------------------------------------\n"	
					arr.append([])
					set_param(DEVICE, SHARED_MEMORY, DOUBLE_PRECISION, MATRIX_LOW, MATRIX_STEP, MATRIX_UP, AVERAGE_RUNS, LU)
					os.system("../step-1-build-desktop/step-1")
					result = open(Filename(DEVICE,SHARED_MEMORY,DOUBLE_PRECISION,LU), 'r')
					result_list = result.readlines()
					result.close()
					
					for line in result_list:
						if(str.find(line,"n rows CPU factorization") == -1):
							numbers = map(float, line.split())
							arr[run].append(numbers)
					if run==0:
						global_arr=[]
						global_std=[]
						for i in range(0,len(arr[run])):
							global_arr.append([])	
							global_std.append([])				
							for j in range(0,len(arr[run][i])):
								global_arr[i].append(0)
								global_std[i].append(0)
	
					for i in range(0,len(arr[run])):
						for j in range(0,len(arr[run][i])):
							global_arr[i][j]+=arr[run][i][j]
						
				#	print "Run \n"
				#	print arr[run]
				for i in range(0,len(arr[run])):
					for j in range(0,len(arr[run][i])):
						global_arr[i][j]*=1.0/AVERAGE_RUNS
	
				#print "\n Averaged data Matrix"			
				#print global_arr
				#print global_std
				#print "\n"
			
	
				for run in range(0,AVERAGE_RUNS):						
					for i in range(0,len(arr[run])):
						for j in range(0,len(arr[run][i])):
							if run==0:
								global_std[i][j]=0							
							#	a=0
							global_std[i][j]+=(arr[run][i][j]-global_arr[i][j])*(arr[run][i][j]-global_arr[i][j])
				print global_arr
				print "\n"
				print global_std
				print "\n"
				#print arr[0]
				print "\n"
				for i in range(0,len(arr[run])):
					for j in range(0,len(arr[run][i])):		
						global_std[i][j]=(global_std[i][j]/AVERAGE_RUNS)**0.5							
				#global_arr=global_arr+global_std	
	
				#print global_arr
				#print "\n"
	
				
				gnu_arr = []
				for i in range(0, len(global_arr)):
					gnu_arr.append([])
					for j in range(0, len(global_arr[i])):
						gnu_arr[i].append(global_arr[i][j])
						gnu_arr[i].append(global_std[i][j]/(AVERAGE_RUNS-1))			
				
	
				averaged = open(Filename(DEVICE,SHARED_MEMORY,DOUBLE_PRECISION,LU)+"_avg", 'w')
				for item in gnu_arr:
					for value in item:
	  					averaged.write("%.3f\t" % value)
					averaged.write("\n")
	

	
				#averaged.close()
							
				
					

