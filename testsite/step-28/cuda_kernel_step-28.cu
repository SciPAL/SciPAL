// To outfox QTCreator's syntax highlighting and especially
// nvcc we put all cuda-related code into files with names
// ending on .cu.c and include them here.
// Then, in the project file we only have to take of one source
// file. This reduces the amount of maintenance.

#include "cuda_utils.cu.c"
#include "cuda_kernel_step-28.cu.c"
