# Gathering *.cu and *.cpp files into one monolithic *.cc file so that 
# the deal.II doc generator can be used to create a documentation for a cuda project
# author : Stephan Kramer, 11.12.2009, Institut fuer Numerische und angewandte Mathematik, U Goettingen 

from glob import iglob
import shutil
import os
import sys


def copy_from_directory_content(target_dir, dst_filename, destination, tmp_br_filename):
# get all cuda-specific files
	kernels = []
	cu_file = []

	for filename in iglob(os.path.join('.', target_dir, '*.cu*')):
		result = os.path.splitext(filename)
		print "found cu-file : " + str(result[0])
		if result[0].endswith("kernel.cu"):
			kernels.append(filename)
		else:
			cu_file.append(filename)
	

#  dump kernels
	
	

	for filename in kernels:
		print "gathering cu kernel " + filename
		shutil.copyfileobj(open(filename, 'rb'), destination)
		shutil.copyfileobj(open(tmp_br_filename, 'rb' ), destination)

	if (len(cu_file) > 1000):
		print "Too many cu-files"
		os.remove(dst_filename)
		os.remove(tmp_br_filename)
		print "abort."
		return

	for filename in cu_file:
		print "gathering cu driver " + filename
		shutil.copyfileobj(open(filename, 'rb'), destination)
		shutil.copyfileobj(open(tmp_br_filename, 'rb' ), destination)

# dump C++ sources
	for filename in iglob(os.path.join(target_dir, '*.h')):
		print "gathering header " + filename
		shutil.copyfileobj(open(filename, 'rb'), destination)
		shutil.copyfileobj(open(tmp_br_filename, 'rb' ), destination)
		
	for filename in iglob(os.path.join(target_dir, '*.cpp')):
		print "gathering source " + filename
		shutil.copyfileobj(open(filename, 'rb'), destination)
		shutil.copyfileobj(open(tmp_br_filename, 'rb' ), destination)


def copy_from_file_list(target_dir, destination, tmp_br_filename):
	oldcwd = os.getcwd()
	os.chdir(target_dir)
	print "copy_from_file_list cwd : " + os.getcwd()
	for filename in open(os.path.join(target_dir, 'doxygen_filelist')):
		if not len(filename.strip('\n').strip(' ')) > 0:
			continue
		if filename.startswith('#'):
			continue
		shutil.copyfileobj(open(filename.strip('\n').strip(' '), 'rb'), destination)
		shutil.copyfileobj(open(tmp_br_filename, 'rb' ), destination)
	os.chdir(oldcwd)


def cat_cu_files():
	target_dir = sys.argv[1]
	dst_filename = os.path.join(target_dir, 'step-cu.cc')
	destination = open(dst_filename,'wb')

	# dummy file for storing some line breaks
	tmp_br_filename = os.path.join('.', target_dir, 'my_br.txt')
	tmp_br = open(tmp_br_filename, 'w' )
	tmp_br.write("\n\n\n")
	tmp_br.close()

	if os.path.exists(os.path.join(target_dir, 'doxygen_filelist')):
		copy_from_file_list(target_dir, destination, tmp_br_filename)
	else:
		copy_from_directory_content(target_dir, dst_filename, destination, tmp_br_filename)

	destination.close()
	os.remove(tmp_br_filename)


cat_cu_files()
