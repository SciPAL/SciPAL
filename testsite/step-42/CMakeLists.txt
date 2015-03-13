include_directories(. ..)

cuda_add_executable(
  step-42
  cuda_kernel_step-42.cu
  SimParams.cpp
  step-42.cpp
)
CUDA_ADD_CUBLAS_TO_TARGET(step-42)

target_link_libraries(step-42 Qt4::QtCore deal_II ${OPENBLAS_LIBRARY})
