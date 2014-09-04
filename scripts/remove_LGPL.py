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

#!/usr/bin/python
# usage: python remove_LGPL.py DICT
# where DICT is the path relative to this script. the LGPL header is removed from all files in that folder.
import os
import string
import sys
from time import *


path = sys.argv[1]
print str(path)

dirList=os.listdir(path)

#get number of lines
lines = 0
for line in open("file_header.txt"):
    lines += 1

print lines


for fname in dirList:
    fullpath = path + fname
    os.system("tail -n +" + str(lines) + " " + str(fullpath) +  " > temp.txt")
    os.system("mv  temp.txt " + str(fullpath))
    print "tail -n +" + str(lines) + " " + str(fullpath) +  " > temp.txt"





