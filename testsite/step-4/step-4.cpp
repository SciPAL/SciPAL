
#include "step-4.hh"

// @sect4{Function: main}
//
// Here the user input is processed and hand overed  to the tests.
int main(int argc, char *argv[])
{
    cudaGetDeviceCount(&global_data.n_CUDA_devices);
    std::cout
            << "N available CUDA devices : "
            << global_data.n_CUDA_devices
            << std::endl;

    step4::run_qr_tests(argc, argv);

    std::cout << std::endl;
}
