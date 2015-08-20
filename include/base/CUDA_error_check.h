#ifndef CUDA_ERROR_CHECK_H
#define CUDA_ERROR_CHECK_H

#include <QtGlobal>

#include <stdio.h>
#include <execinfo.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cuda_runtime.h>
#ifndef Q_OS_OSX
#include <sys/prctl.h>
#endif
//! code from
//! http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__, __FUNCTION__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, const char *func,bool abort=true)
{
#ifndef DEBUG
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s line: %d\n function:%s\n", cudaGetErrorString(code), file, line, func);

#ifndef Q_OS_OSX
      //print backtrace
      char pid_buf[30];
      sprintf(pid_buf, "%d", getpid());
      char name_buf[512];
      name_buf[readlink("/proc/self/exe", name_buf, 511)]=0;
      prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY, 0, 0, 0);
      int child_pid = fork();
      if (!child_pid) {
          dup2(2,1); // redirect output to stderr
          fprintf(stdout,"stack trace for %s pid=%s\n",name_buf,pid_buf);
          execlp("gdb", "gdb", "--batch", "-n", "-ex", "thread", "-ex", "bt", name_buf, pid_buf, NULL);
      } else {
          waitpid(child_pid,NULL,0);
      }
#endif//Q_OS_OSX

      if (abort) exit(code);
   }
#endif//DEBUG
}

#endif // CUDA_ERROR_CHECK_H
