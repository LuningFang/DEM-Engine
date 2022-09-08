// DEM integration related custom kernels
#include <kernel/DEMHelperKernels.cu>
#include <DEM/DEMDefines.h>

// Apply presecibed velocity and report whether the `true' physics should be skipped, rather than added on top of that
template <typename T1, typename T2>
inline __device__ void applyPrescribedVel(bool& LinPrescribed,
                                          bool& RotPrescribed,
                                          T1& vX,
                                          T1& vY,
                                          T1& vZ,
                                          T2& omgBarX,
                                          T2& omgBarY,
                                          T2& omgBarZ,
                                          const smug::family_t& family,
                                          const float& t) {
    switch (family) {
        _velPrescriptionStrategy_;
        default:
            LinPrescribed = false;
            RotPrescribed = false;
    }
}

// Apply presecibed location and report whether the `true' physics should be skipped, rather than added on top of that
template <typename T1, typename T2>
inline __device__ void applyPrescribedPos(bool& LinPrescribed,
                                          bool& RotPrescribed,
                                          T1& X,
                                          T1& Y,
                                          T1& Z,
                                          T2& oriQw,
                                          T2& oriQx,
                                          T2& oriQy,
                                          T2& oriQz,
                                          const smug::family_t& family,
                                          const float& t) {
    switch (family) {
        _posPrescriptionStrategy_;
        default:
            LinPrescribed = false;
            RotPrescribed = false;
    }
}

inline __device__ void integrateVel(smug::bodyID_t thisClump,
                                    smug::DEMSimParams* simParams,
                                    smug::DEMDataDT* granData,
                                    float h,
                                    float t) {
    smug::family_t family_code = granData->familyID[thisClump];
    bool LinPrescribed = false, RotPrescribed = false;
    applyPrescribedVel<float, float>(LinPrescribed, RotPrescribed, granData->vX[thisClump], granData->vY[thisClump],
                                     granData->vZ[thisClump], granData->omgBarX[thisClump],
                                     granData->omgBarY[thisClump], granData->omgBarZ[thisClump], family_code, (float)t);

    if (!LinPrescribed) {
        granData->vX[thisClump] += (granData->aX[thisClump] + simParams->Gx) * h;
        granData->vY[thisClump] += (granData->aY[thisClump] + simParams->Gy) * h;
        granData->vZ[thisClump] += (granData->aZ[thisClump] + simParams->Gz) * h;
    }

    if (!RotPrescribed) {
        granData->omgBarX[thisClump] += granData->alphaX[thisClump] * h;
        granData->omgBarY[thisClump] += granData->alphaY[thisClump] * h;
        granData->omgBarZ[thisClump] += granData->alphaZ[thisClump] * h;
    }
}

inline __device__ void locateNewVoxel(smug::voxelID_t& voxel, int64_t& locX_tmp, int64_t& locY_tmp, int64_t& locZ_tmp) {
    smug::voxelID_t voxelX;
    smug::voxelID_t voxelY;
    smug::voxelID_t voxelZ;
    IDChopper<smug::voxelID_t, smug::voxelID_t>(voxelX, voxelY, voxelZ, voxel, _nvXp2_, _nvYp2_);

    // DEM_MAX_SUBVOXEL is int64 and large enough to handle DEM_VOXEL_RES_POWER2 == 16 or 32
    voxelX += div_floor<int64_t, int64_t>(locX_tmp, smug::DEM_MAX_SUBVOXEL);
    voxelY += div_floor<int64_t, int64_t>(locY_tmp, smug::DEM_MAX_SUBVOXEL);
    voxelZ += div_floor<int64_t, int64_t>(locZ_tmp, smug::DEM_MAX_SUBVOXEL);
    locX_tmp = mod_floor<int64_t, int64_t>(locX_tmp, smug::DEM_MAX_SUBVOXEL);
    locY_tmp = mod_floor<int64_t, int64_t>(locY_tmp, smug::DEM_MAX_SUBVOXEL);
    locZ_tmp = mod_floor<int64_t, int64_t>(locZ_tmp, smug::DEM_MAX_SUBVOXEL);

    // TODO: Should add a check here where, if negative voxel component spotted, stop the simulation

    IDPacker<smug::voxelID_t, smug::voxelID_t>(voxel, voxelX, voxelY, voxelZ, _nvXp2_, _nvYp2_);
}

inline __device__ void integratePos(smug::bodyID_t thisClump, smug::DEMDataDT* granData, float h, float t) {
    // Location accuracy is up to integer level anyway
    int64_t locX_tmp = (int64_t)granData->locX[thisClump];
    int64_t locY_tmp = (int64_t)granData->locY[thisClump];
    int64_t locZ_tmp = (int64_t)granData->locZ[thisClump];

    smug::family_t family_code = granData->familyID[thisClump];
    bool LinPrescribed = false, RotPrescribed = false;
    applyPrescribedPos<int64_t, smug::oriQ_t>(
        LinPrescribed, RotPrescribed, locX_tmp, locY_tmp, locZ_tmp, granData->oriQw[thisClump],
        granData->oriQx[thisClump], granData->oriQy[thisClump], granData->oriQz[thisClump], family_code, (float)t);

    if (!LinPrescribed) {
        locX_tmp += (int64_t)((double)granData->vX[thisClump] / _l_ * h);
        locY_tmp += (int64_t)((double)granData->vY[thisClump] / _l_ * h);
        locZ_tmp += (int64_t)((double)granData->vZ[thisClump] / _l_ * h);
        smug::voxelID_t newVoxel = granData->voxelID[thisClump];
        locateNewVoxel(newVoxel, locX_tmp, locY_tmp, locZ_tmp);
        granData->voxelID[thisClump] = newVoxel;
        granData->locX[thisClump] = locX_tmp;
        granData->locY[thisClump] = locY_tmp;
        granData->locZ[thisClump] = locZ_tmp;
    }

    if (!RotPrescribed) {
        // Then integrate the quaternion
        // Exp map-based rotation angle calculation
        double deltaQ0 = 1;
        double deltaQ1 = granData->omgBarX[thisClump];
        double deltaQ2 = granData->omgBarY[thisClump];
        double deltaQ3 = granData->omgBarZ[thisClump];
        double len = sqrt(deltaQ1 * deltaQ1 + deltaQ2 * deltaQ2 + deltaQ3 * deltaQ3);
        double theta = 0.5 * h * len;  // 0.5*dt*len, delta rotation
        if (len > 0) {
            deltaQ0 = cos(theta);
            double s = sin(theta) / len;
            deltaQ1 *= s;
            deltaQ2 *= s;
            deltaQ3 *= s;
        }
        // Note: Yes it is Quat * deltaRot, not the other way around. Also, Hamilton product should automatically
        // maintain the unit-ness of quaternions.
        HamiltonProduct<float>(granData->oriQw[thisClump], granData->oriQx[thisClump], granData->oriQy[thisClump],
                               granData->oriQz[thisClump], granData->oriQw[thisClump], granData->oriQx[thisClump],
                               granData->oriQy[thisClump], granData->oriQz[thisClump], deltaQ0, deltaQ1, deltaQ2,
                               deltaQ3);
    }
}

__global__ void integrateClumps(smug::DEMSimParams* simParams, smug::DEMDataDT* granData, float t) {
    smug::bodyID_t thisClump = blockIdx.x * blockDim.x + threadIdx.x;
    if (thisClump < simParams->nOwnerBodies) {
        integrateVel(thisClump, simParams, granData, simParams->h, t);
        integratePos(thisClump, granData, simParams->h, t);
    }
}
