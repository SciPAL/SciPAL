#!/bin/bash

make -f Makefile.dox_doc clean
cd ../../doc/autogen/CUDA_HPC_Praktikum
rm -rf *
cd ../../tutorial
rm -rf doxygen/* generated/*

