#ifndef CUDA_HELPERS_CUH
#define CUDA_HELPERS_CUH
#include <cuda.h>
#include <cuda_runtime_api.h>

// Macro to catch CUDA errors in CUDA runtime calls
#define CUDA_SAFE_CALL(ans)                   \
    {                                         \
        gpuAssert((ans), __FILE__, __LINE__); \
    }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s in %s, line %d\n", cudaGetErrorString(code), file, line);
        if (abort)
            exit(code);
    }
}
#endif