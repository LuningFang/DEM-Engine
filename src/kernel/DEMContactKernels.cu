// DEM contact detection related custom kernels
#include <granular/DataStructs.h>
#include <granular/GranularDefines.h>

__global__ void getNumberOfBinsEachSphereTouches(sgps::DEMSimParams* simParams, sgps::DEMDataKT* granData) {
    // __shared__ const distinctSphereRadii[???];
    sgps::bodyID_t sphereID = blockIdx.x * blockDim.x + threadIdx.x;
    if (sphereID < simParams->nSpheresGM) {
    }
}

__global__ void populateBinSphereTouchingPairs(sgps::DEMSimParams* simParams, sgps::DEMDataKT* granData) {
    // __shared__ const distinctSphereRadii[???];
    sgps::bodyID_t sphereID = blockIdx.x * blockDim.x + threadIdx.x;
    if (sphereID < simParams->nSpheresGM) {
    }
}
