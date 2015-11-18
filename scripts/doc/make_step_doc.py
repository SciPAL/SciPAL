#!/usr/bin/env python
# Creation of a new step-xy project so that
# the deal.II doc generator can be used to create a documentation for a cuda project
# author : Stephan Kramer, 11.12.2009, Institut fuer Numerische und angewandte Mathematik, U Goettingen 

#from glob import iglob
from glob import glob
import shutil
import os
import sys



def prepare_makefile():

    if (len(sys.argv) < 2):
        print "Usage : Provide the numbers of steps for which a documentation\nshould be generated as colon separated list of ranges and numbers.\ne.g.: " + str(sys.argv[0]) + " 3-6,8,11"
        return

    raw_step_numbers = sys.argv[1:]
    print raw_step_numbers

    tmp_step_numbers = []
    for e in raw_step_numbers:
	aux = e.split(",")
	for a in aux:
	    print a
	    if a.find("-") > -1:
		[a_min, a_max] = a.split("-")
		if int(a_min) > int(a_max):
		    tmp = a_max
		    a_max = a_min
		    a_min = tmp
		for s in range(int(a_min),int(a_max)+1):
		    tmp_step_numbers.append(s)
	    else:
		step_no = int(0)
		try:
		    step_no = int(a)
		except:
		    step_no = 0
		if (step_no > 0):
		    tmp_step_numbers.append(step_no)

    print tmp_step_numbers
    tmp_step_numbers.sort()
    
    step_numbers = []
	# remove duplicates
    for n in tmp_step_numbers:
	if not n in step_numbers:
	    step_numbers.append(n) 
    
    print tmp_step_numbers
    print step_numbers
    
    makefile = open("../../doc/tutorial/Makefile",'r')
    lines = makefile.readlines()
    makefile.close()

                # Ersetze namespace-Namen, etc ...
    makefile = open("../../doc/tutorial/Makefile",'w')
    for line in lines:
        if line.startswith("example-names"):
	    tmpl = line.split("=")
	    aux = tmpl[0] + " = "
	    steps = ""
	    for s in step_numbers:
		steps = steps + "step-%d " % s
            line = aux + steps + "\n"
#	    print line
        makefile.write(line)
    makefile.close()    
    return


def create_step_doc():
    os.system("make -f Makefile.dox_doc clean")
    os.system("make -f Makefile.dox_doc")


def repair_bogus_step_names():
    current_wd = os.getcwd()
    os.chdir("../../doc/autogen/CUDA_HPC_Praktikum/")
    step_names = glob("step__*.html")
    for n in step_names:
	number = n.split(".")[0].split("__")[1]
	print "bogus step name : " + number
        src = n
        dst = "step_"+number+".html"
        os.rename(src, dst)
    os.chdir(current_wd)

# main
prepare_makefile()
create_step_doc()
repair_bogus_step_names()


