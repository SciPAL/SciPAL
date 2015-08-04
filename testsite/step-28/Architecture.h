// @sect3{File: Architecture.h}
#ifndef STEP28ARCHITECTURE_H
#define STEP28ARCHITECTURE_H

// @sect4{enum: Architecture}
//
// Enum to choose the architecure at runtime. @p both is for debugging: the values
// calculated by CUDA are checked against those calculated on the CPU
namespace step28 {

enum Architecture { cpu, cuda, both };

}

#endif // PARALLELARCH_H
