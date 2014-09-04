import os
from time import sleep

var = raw_input("Which setting shall be tested? '(d)evice', '(m)emory', '(p)recision (between double and float)' or 'decomposition (t)ype (between cholesky and LU)' ")

if (var.upper() == "D"):
	setting = "device"
elif (var.upper() == "M"):
	setting = "memory"
elif (var.upper() == "P"):	
	setting = "prec"
elif (var.upper() == "T"):
	setting = "decompType"

average_runs = raw_input("How many averaging runs shall be taken? ")

print("Wait for simulations to start...\n");
sleep(2)
os.system("python multi_tests3.py " + average_runs )
os.system("python plotResults2.py "+setting)
os.system("okular results.dvi")
 
