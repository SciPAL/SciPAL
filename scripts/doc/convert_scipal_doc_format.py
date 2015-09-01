#!/usr/bin/env python
# Creation of html doc pages for SciPal library
# author : Stephan Kramer, 25.05.2012, Institut fuer Numerische und angewandte Mathematik, U Goettingen,
# J. Hagemann, same affiliation 


#from glob import iglob
import glob
import shutil
import os
import sys
import string




#get all scipal files
filenames = glob.glob("../../include/lac/*.h*")
filenames = filenames + glob.glob("../../include/base/*.h*")
# ...
print filenames

# // -> //!
for fname in filenames:
	fl = open(fname,'r')
	lines = fl.readlines()
	fl.close()
	fl = open(fname,'w')
	for line in lines:
		if(str.find(line,"//") != -1 and str.find(line,"@sect") == -1):
			tmp = line.replace("//", "//!")
			line = tmp
		fl.write(line)
#print line
#	fl = open(fname,'w')
#	fl.writelines(lines)
	fl.close()

