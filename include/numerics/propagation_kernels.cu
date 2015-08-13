//! \file
// Um QTCreators mangelndes Syntaxhighlighting fuer cu-files
// auszutricksen werden die eigentlichen Treiber-Klassen in eine 
// C-header-Datei geschrieben und diese hier eingebunden, so dass
// nvcc denkt, er haette ein cu-file for sich.

#include "propagation_kernels.cu.c"

