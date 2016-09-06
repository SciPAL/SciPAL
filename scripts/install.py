#!/usr/bin/python
#This file is part of SciPAL.

#    SciPAL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    SciPAL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public License
#    along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.

#Copyright  S. C. Kramer , J. Hagemann  2010 - 2014

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import string
import sys
from time import *

print ("This is the SciPAL install script \n"
        "We need a few paths from you for the initial configuration. \n"
        "Problems? : scipal@googlegroups.com");

#get the location of install.py
script_path = os.path.dirname(os.path.realpath(__file__));
#change in scripts directory, so that the relative paths work
os.chdir(script_path);

dealii = raw_input("Enter path to deal.ii:");


conf = open("../config/dealii_simple_conf.pro_in","r")

conf_list = conf.readlines()
conf.close()
i = 0
for line in conf_list:
    if(string.find(line,"DEALHOME =") != -1):
        conf_list[i] = "DEALHOME ="+str(dealii)+"\n";
    i +=1;

conf = open("../config/dealii_simple_conf.pro","w")
conf.writelines(conf_list)
conf.close()

#CUDA config
correct = False
while not correct:
	print "Is your CUDA installation located here: /usr/local/cuda?";
	answer = raw_input("(y)es/(n)")
	if answer:
		if answer[0].lower() == "y" :
			print(answer[0].lower());
			dir = "/usr/local/cuda"
			correct = True;

		if answer[0].lower() == "n":
    			answer = raw_input("Enter path to CUDA installation:");
			dir = answer
			correct = True;


conf = open("../config/scipal_conf.pro_in","r")
conf_list = conf.readlines()
conf.close()
i = 0
for line in conf_list:
    if(str.find(line,"CUDA_DIR = /usr/local/cuda") != -1):
        conf_list[i] = "CUDA_DIR = "+str(dir)+"\n";

    i +=1;

conf = open("../config/scipal_conf.pro","w")
conf.writelines(conf_list)
conf.close()

conf = open("doc/options.136","r")
conf_list = conf.readlines()
conf.close()

i = 0
for line in conf_list:
    if(string.find(line,"INCLUDE_PATH           =") != -1):
        conf_list[i] = "INCLUDE_PATH           ="+str(dealii)+"/base "+str(dealii)+"/lac "+str(dealii)+"/ "+str(dir)+"/include\n";
    i +=1;

conf = open("doc/options.136","w")
conf.writelines(conf_list)
conf.close()



##further deal ii config
#petsc
#mpi
#metis
#tbb

print ("You can rerun this script to change paths. \n");


