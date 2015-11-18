#!/bin/bash

make -f Makefile.dox_doc clean
./make_step_doc.py 1,2,4,11-18
make -f Makefile.dox_doc
