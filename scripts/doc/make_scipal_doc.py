#!/usr/bin/env python
# Creation of html doc pages for SciPal library
# author : Stephan Kramer, 25.05.2012, Institut fuer Numerische und angewandte Mathematik, U Goettingen,
# J. Hagemann, same affiliation 


import shutil
import os
import sys
import string

if os.path.exists("../../doc/autogen/SciPAL"):
	shutil.rmtree("../../doc/autogen/SciPAL")
os.mkdir("../../doc/autogen/SciPAL")

#call doxygen
os.system("doxygen scipal.dox")




