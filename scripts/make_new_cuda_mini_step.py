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

#Copyright  S. C. Kramer , J. Hagemann, T. Koehler  2010 - 2015


# Creation of a new step-xy mini project so that
# the deal.II doc generator can be used to create a documentation for a cuda project
# author : Stephan Kramer, 11.12.2009, Institut fuer Numerische und angewandte Mathematik, U Goettingen 

#from glob import iglob
from glob import glob
import shutil
import os
import sys



def copy_template():

    if (len(sys.argv) != 2):
        print "Usage : " + str(sys.argv[0]) + " <number>-<number>"
        print "if only one <number> is given it indicates the number of the step"
        print "<number>-<number> generates a contiguous set of new steps"
        return

    step_numbers = sys.argv[1].split("-")
    step_min = int(0)
    step_max = int(0)
    if len(step_numbers) > 2:
        print "illegal input"
        return
    elif len(step_numbers) == 2:
        step_min = int(step_numbers[0])
        step_max = int(step_numbers[1])+1
    elif len(step_numbers) == 1:
        step_min = int(step_numbers[0])
        step_max = step_min+1
    print "step range : %d - %d " % (step_min, step_max) 
    
    print "cwd : " + os.getcwd() 	
    while str(os.path.split(os.getcwd())[1]) != "SciPAL":
        os.chdir("..")
        print "cwd : " + os.getcwd()
    os.chdir("testsite")
    print "cwd : " + os.getcwd() 	
    
    for step in range(step_min,step_max):
        step_dirname = "mini-step-%d" % int(step)
        print "new step-dir : " + step_dirname
        print "Start cwd : " + os.getcwd()
        if os.path.exists(step_dirname):
            print "mini-step %d already exists" % step
        else:
            shutil.copytree("step-mini", step_dirname)
            os.chdir(step_dirname)
            print "Start : " + os.getcwd()
            template_files = glob("*mini*.*")
            template_files.append(glob("doxygen_filelist")[0])
            for f in template_files:
                file = open(f,'r')
                print "file to process : %s" % str(f)
                lines = file.readlines()
                file.close()
                # Ersetze namespace-Namen, etc ...
                file = open(f,'w')
                for line in lines:
                    if line.find("stepmini") > 1:
                        tmpl = line.replace("stepmini", "ministep%d" % step)
                        line = tmpl
                    if line.startswith("stepmini"):
                        tmpl = line.replace("stepmini", "ministep%d" % step)
                        line = tmpl
                    if line.find("MINI") > 1:
                        tmpl = line.replace("MINI", "MINI_%d" % step)
                        line = tmpl
                    if line.find("step-mini") > 1:
                        tmpl = line.replace("step-mini", "mini-step-%d" % step)
                        line = tmpl
                    if line.startswith("step-mini"):
                        tmpl = line.replace("step-mini", "mini-step-%d" % step)
                        line = tmpl
                    if line.find("step_mini") > 1:
                        tmpl = line.replace("step_mini", "mini_step_%d" % step)
                        line = tmpl
                    file.write(line)
                file.close()
                if f != "doxygen_filelist":
                    shutil.move(f, f.replace("step-mini", "mini-step-%d" % step))
            # clean up
            print "Done : " + os.getcwd()
            if os.path.exists("mini-step-" + str(step) + ".pro.user"):
                os.remove("mini-step-" + str(step) + ".pro.user")
            if os.path.exists("autoInstantiations.h"):
                os.remove("autoInstantiations.h")

            shutil.rmtree(".svn")
            shutil.rmtree("prm/.svn")
            shutil.rmtree("doc/.svn")
            shutil.rmtree("doc/images/.svn")
            os.chdir("..")


copy_template()
