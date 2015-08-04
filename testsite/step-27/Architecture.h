// @sect3{File: Architecture.h}
#ifndef STEP27ARCHITECTURE_H
#define STEP27ARCHITECTURE_H

// @sect4{enum: Architecture}
//
// Enum to choose the architecure at runtime. @p both is for debugging: the values
// calculated by CUDA are checked against those calculated on the CPU
namespace step27 {

enum Architecture { cpu, cuda, both };

}

#endif // PARALLELARCH_H
