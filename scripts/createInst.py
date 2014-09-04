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

import sys, time, re, datetime, os

MakefileVersion = sys.argv[1]

lines=sys.stdin.readlines()

newImpl = []    # Save all new instantiations
newError = []   # Save all other errors

for line in lines:  # Go through the whole compiler output
        matchobject = re.search("\`void Impl(.*)\'", line)  # Look for "`void Impl( * '"
	if (matchobject!=None):
                newImpl.append("template void Impl"+matchobject.groups(1)[0]+";")   # if it was found a new instantiation is needed
        else:
                newError.append(line)   # otherwise it is another error
		#print " "
if (len(newImpl)!=0):   # if a new instantiation is needed print them into stdout
	dt = datetime.datetime.now()
	print "// The following "+str(len(set([x for x in newImpl])))+" instantiations are automatically generated ["+dt.isoformat(' ')+"]:"
	for line in set([x for x in newImpl]):
		print line
        sys.stdout.flush()  # important, otherwise the buffer is not empty at the right time
        os.system("make clean > /dev/null");    # very important: otherwise the include is not redone!
        os.system("make -f Makefile."+MakefileVersion+" >&2");  # second build, now with the new instantiations.
else:
        for line in newError:   # print all errors if there were no new instantiations are created.
                sys.stderr.write(line)
        sys.stderr.write("\n")
