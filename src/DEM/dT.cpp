//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <cstring>
#include <iostream>
#include <thread>
#include <algorithm>

#include <chpf.hpp>
#include <core/ApiVersion.h>
#include <core/utils/JitHelper.h>
#include <DEM/dT.h>
#include <DEM/kT.h>
#include <DEM/HostSideHelpers.hpp>
#include <nvmath/helper_math.cuh>
#include <DEM/DEMDefines.h>

#include <algorithms/DEMCubBasedSubroutines.h>

namespace sgps {

// Put sim data array pointers in place
void DEMDynamicThread::packDataPointers() {
    granData->inertiaPropOffsets = inertiaPropOffsets.data();
    granData->familyID = familyID.data();
    granData->voxelID = voxelID.data();
    granData->locX = locX.data();
    granData->locY = locY.data();
    granData->locZ = locZ.data();
    granData->aX = aX.data();
    granData->aY = aY.data();
    granData->aZ = aZ.data();
    granData->vX = vX.data();
    granData->vY = vY.data();
    granData->vZ = vZ.data();
    granData->oriQ0 = oriQ0.data();
    granData->oriQ1 = oriQ1.data();
    granData->oriQ2 = oriQ2.data();
    granData->oriQ3 = oriQ3.data();
    granData->omgBarX = omgBarX.data();
    granData->omgBarY = omgBarY.data();
    granData->omgBarZ = omgBarZ.data();
    granData->alphaX = alphaX.data();
    granData->alphaY = alphaY.data();
    granData->alphaZ = alphaZ.data();
    granData->idGeometryA = idGeometryA.data();
    granData->idGeometryB = idGeometryB.data();
    granData->contactType = contactType.data();
    granData->familyMasks = familyMaskMatrix.data();

    // granData->idGeometryA_buffer = idGeometryA_buffer.data();
    // granData->idGeometryB_buffer = idGeometryB_buffer.data();
    // granData->contactType_buffer = contactType_buffer.data();
    // granData->contactMapping_buffer = contactMapping_buffer.data();

    granData->contactForces = contactForces.data();
    granData->contactTorque_convToForce = contactTorque_convToForce.data();
    granData->contactPointGeometryA = contactPointGeometryA.data();
    granData->contactPointGeometryB = contactPointGeometryB.data();
    granData->contactHistory = contactHistory.data();
    granData->contactDuration = contactDuration.data();

    // The offset info that indexes into the template arrays
    granData->ownerClumpBody = ownerClumpBody.data();
    granData->clumpComponentOffset = clumpComponentOffset.data();
    granData->clumpComponentOffsetExt = clumpComponentOffsetExt.data();
    granData->materialTupleOffset = materialTupleOffset.data();

    // Template array pointers
    granData->radiiSphere = radiiSphere.data();
    granData->relPosSphereX = relPosSphereX.data();
    granData->relPosSphereY = relPosSphereY.data();
    granData->relPosSphereZ = relPosSphereZ.data();
    granData->massOwnerBody = massOwnerBody.data();
    granData->mmiXX = mmiXX.data();
    granData->mmiYY = mmiYY.data();
    granData->mmiZZ = mmiZZ.data();
    granData->EProxy = EProxy.data();
    granData->nuProxy = nuProxy.data();
    granData->CoRProxy = CoRProxy.data();
    granData->muProxy = muProxy.data();
    granData->CrrProxy = CrrProxy.data();
}

void DEMDynamicThread::packTransferPointers(DEMKinematicThread*& kT) {
    // These are the pointers for sending data to dT
    granData->pKTOwnedBuffer_voxelID = kT->granData->voxelID_buffer;
    granData->pKTOwnedBuffer_locX = kT->granData->locX_buffer;
    granData->pKTOwnedBuffer_locY = kT->granData->locY_buffer;
    granData->pKTOwnedBuffer_locZ = kT->granData->locZ_buffer;
    granData->pKTOwnedBuffer_oriQ0 = kT->granData->oriQ0_buffer;
    granData->pKTOwnedBuffer_oriQ1 = kT->granData->oriQ1_buffer;
    granData->pKTOwnedBuffer_oriQ2 = kT->granData->oriQ2_buffer;
    granData->pKTOwnedBuffer_oriQ3 = kT->granData->oriQ3_buffer;
    granData->pKTOwnedBuffer_familyID = kT->granData->familyID_buffer;
}

void DEMDynamicThread::changeFamily(unsigned int ID_from, unsigned int ID_to) {
    family_t ID_from_impl = ID_from;
    family_t ID_to_impl = ID_to;
    std::replace_if(
        familyID.begin(), familyID.end(), [ID_from_impl](family_t& i) { return i == ID_from_impl; }, ID_to_impl);
}

void DEMDynamicThread::setSimParams(unsigned char nvXp2,
                                    unsigned char nvYp2,
                                    unsigned char nvZp2,
                                    float l,
                                    double voxelSize,
                                    double binSize,
                                    binID_t nbX,
                                    binID_t nbY,
                                    binID_t nbZ,
                                    float3 LBFPoint,
                                    float3 G,
                                    double ts_size,
                                    float expand_factor,
                                    float approx_max_vel,
                                    float expand_safety_param) {
    simParams->nvXp2 = nvXp2;
    simParams->nvYp2 = nvYp2;
    simParams->nvZp2 = nvZp2;
    simParams->l = l;
    simParams->voxelSize = voxelSize;
    simParams->binSize = binSize;
    simParams->LBFX = LBFPoint.x;
    simParams->LBFY = LBFPoint.y;
    simParams->LBFZ = LBFPoint.z;
    simParams->Gx = G.x;
    simParams->Gy = G.y;
    simParams->Gz = G.z;
    simParams->h = ts_size;
    simParams->beta = expand_factor;
    simParams->approxMaxVel = approx_max_vel;
    simParams->expSafetyParam = expand_safety_param;
    simParams->nbX = nbX;
    simParams->nbY = nbY;
    simParams->nbZ = nbZ;
}

float DEMDynamicThread::getKineticEnergy() {
    // // We can use temp vectors as we please
    // size_t quarryTempSize = (size_t)simParams->nOwnerBodies * sizeof(double);
    // double* KEArr = (double*)stateOfSolver_resources.allocateTempVector(1, quarryTempSize);
    // size_t returnSize = sizeof(double);
    // double* KE = (double*)stateOfSolver_resources.allocateTempVector(2, returnSize);
    // size_t blocks_needed_for_KE =
    //     (simParams->nOwnerBodies + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) / SGPS_DEM_MAX_THREADS_PER_BLOCK;
    // quarry_stats_kernels->kernel("computeKE")
    //     .instantiate()
    //     .configure(dim3(blocks_needed_for_KE), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0, streamInfo.stream)
    //     .launch(granData, simParams->nOwnerBodies, KEArr);
    // GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
    // // displayArray<double>(KEArr, simParams->nOwnerBodies);
    // doubleSumReduce(KEArr, KE, simParams->nOwnerBodies, streamInfo.stream, stateOfSolver_resources);
    // return *KE;
    return 0;
}

void DEMDynamicThread::changeOwnerSizes(const std::vector<bodyID_t>& IDs, const std::vector<float>& factors) {
    // Set the gpu for this thread
    // cudaSetDevice(streamInfo.device);
    // cudaStream_t new_stream;
    // cudaStreamCreate(&new_stream);

    // First get IDs and factors to device side
    size_t IDSize = IDs.size() * sizeof(bodyID_t);
    bodyID_t* dIDs = (bodyID_t*)stateOfSolver_resources.allocateTempVector(1, IDSize);
    GPU_CALL(cudaMemcpy(dIDs, IDs.data(), IDSize, cudaMemcpyHostToDevice));
    size_t factorSize = factors.size() * sizeof(float);
    float* dFactors = (float*)stateOfSolver_resources.allocateTempVector(2, factorSize);
    GPU_CALL(cudaMemcpy(dFactors, factors.data(), factorSize, cudaMemcpyHostToDevice));

    size_t idBoolSize = (size_t)simParams->nOwnerBodies * sizeof(notStupidBool_t);
    size_t ownerFactorSize = (size_t)simParams->nOwnerBodies * sizeof(float);
    // Bool table for whether this owner should change
    notStupidBool_t* idBool = (notStupidBool_t*)stateOfSolver_resources.allocateTempVector(3, idBoolSize);
    GPU_CALL(cudaMemset(idBool, 0, idBoolSize));
    float* ownerFactors = (float*)stateOfSolver_resources.allocateTempVector(4, ownerFactorSize);
    size_t blocks_needed_for_marking =
        (IDs.size() + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) / SGPS_DEM_MAX_THREADS_PER_BLOCK;

    // Mark on the bool array those owners that need a change
    misc_kernels->kernel("markOwnerToChange")
        .instantiate()
        .configure(dim3(blocks_needed_for_marking), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0, streamInfo.stream)
        .launch(idBool, ownerFactors, dIDs, dFactors, IDs.size());
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));

    // Change the size of the sphere components in question
    size_t blocks_needed_for_changing =
        (simParams->nSpheresGM + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) / SGPS_DEM_MAX_THREADS_PER_BLOCK;
    misc_kernels->kernel("dTModifyComponents")
        .instantiate()
        .configure(dim3(blocks_needed_for_changing), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0, streamInfo.stream)
        .launch(granData, idBool, ownerFactors, simParams->nSpheresGM);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));

    // cudaStreamDestroy(new_stream);
}

void DEMDynamicThread::allocateManagedArrays(size_t nOwnerBodies,
                                             size_t nOwnerClumps,
                                             unsigned int nExtObj,
                                             size_t nTriEntities,
                                             size_t nSpheresGM,
                                             size_t nTriGM,
                                             unsigned int nAnalGM,
                                             unsigned int nMassProperties,
                                             unsigned int nClumpTopo,
                                             unsigned int nClumpComponents,
                                             unsigned int nJitifiableClumpComponents,
                                             unsigned int nMatTuples) {
    // dT buffer arrays should be on dT and this is to ensure that
    GPU_CALL(cudaSetDevice(streamInfo.device));

    // Sizes of these arrays
    simParams->nSpheresGM = nSpheresGM;
    simParams->nTriGM = nTriGM;
    simParams->nAnalGM = nAnalGM;
    simParams->nOwnerBodies = nOwnerBodies;
    simParams->nOwnerClumps = nOwnerClumps;
    simParams->nExtObj = nExtObj;
    simParams->nTriEntities = nTriEntities;
    simParams->nDistinctMassProperties = nMassProperties;
    simParams->nDistinctClumpBodyTopologies = nClumpTopo;
    simParams->nJitifiableClumpComponents = nJitifiableClumpComponents;
    simParams->nDistinctClumpComponents = nClumpComponents;
    simParams->nMatTuples = nMatTuples;

    // Resize to the number of clumps
    SGPS_DEM_TRACKED_RESIZE(familyID, nOwnerBodies, "familyID", 0);
    SGPS_DEM_TRACKED_RESIZE(voxelID, nOwnerBodies, "voxelID", 0);
    SGPS_DEM_TRACKED_RESIZE(locX, nOwnerBodies, "locX", 0);
    SGPS_DEM_TRACKED_RESIZE(locY, nOwnerBodies, "locY", 0);
    SGPS_DEM_TRACKED_RESIZE(locZ, nOwnerBodies, "locZ", 0);
    SGPS_DEM_TRACKED_RESIZE(oriQ0, nOwnerBodies, "oriQ0", 1);
    SGPS_DEM_TRACKED_RESIZE(oriQ1, nOwnerBodies, "oriQ1", 0);
    SGPS_DEM_TRACKED_RESIZE(oriQ2, nOwnerBodies, "oriQ2", 0);
    SGPS_DEM_TRACKED_RESIZE(oriQ3, nOwnerBodies, "oriQ3", 0);
    SGPS_DEM_TRACKED_RESIZE(vX, nOwnerBodies, "vX", 0);
    SGPS_DEM_TRACKED_RESIZE(vY, nOwnerBodies, "vY", 0);
    SGPS_DEM_TRACKED_RESIZE(vZ, nOwnerBodies, "vZ", 0);
    SGPS_DEM_TRACKED_RESIZE(omgBarX, nOwnerBodies, "omgBarX", 0);
    SGPS_DEM_TRACKED_RESIZE(omgBarY, nOwnerBodies, "omgBarY", 0);
    SGPS_DEM_TRACKED_RESIZE(omgBarZ, nOwnerBodies, "omgBarZ", 0);
    SGPS_DEM_TRACKED_RESIZE(aX, nOwnerBodies, "aX", 0);
    SGPS_DEM_TRACKED_RESIZE(aY, nOwnerBodies, "aY", 0);
    SGPS_DEM_TRACKED_RESIZE(aZ, nOwnerBodies, "aZ", 0);
    SGPS_DEM_TRACKED_RESIZE(alphaX, nOwnerBodies, "alphaX", 0);
    SGPS_DEM_TRACKED_RESIZE(alphaY, nOwnerBodies, "alphaY", 0);
    SGPS_DEM_TRACKED_RESIZE(alphaZ, nOwnerBodies, "alphaZ", 0);

    // Resize the family mask `matrix' (in fact it is flattened)
    SGPS_DEM_TRACKED_RESIZE(familyMaskMatrix, (DEM_NUM_AVAL_FAMILIES - 1) * DEM_NUM_AVAL_FAMILIES / 2,
                            "familyMaskMatrix", DEM_DONT_PREVENT_CONTACT);

    // Resize to the number of geometries
    SGPS_DEM_TRACKED_RESIZE(ownerClumpBody, nSpheresGM + nTriGM, "ownerClumpBody", 0);
    SGPS_DEM_TRACKED_RESIZE(materialTupleOffset, nSpheresGM + nTriGM, "materialTupleOffset", 0);

    // For clump component offset, it's only needed if clump components are jitified
    if (solverFlags.useClumpJitify) {
        SGPS_DEM_TRACKED_RESIZE(clumpComponentOffset, nSpheresGM, "clumpComponentOffset", 0);
        // This extended component offset array can hold offset numbers even for big clumps (whereas
        // clumpComponentOffset is typically uint_8, so it may not). If a sphere's component offset index falls in this
        // range then it is not jitified, and the kernel needs to look for it in the global memory.
        SGPS_DEM_TRACKED_RESIZE(clumpComponentOffsetExt, nSpheresGM, "clumpComponentOffsetExt", 0);
        SGPS_DEM_TRACKED_RESIZE(radiiSphere, nClumpComponents, "radiiSphere", 0);
        SGPS_DEM_TRACKED_RESIZE(relPosSphereX, nClumpComponents, "relPosSphereX", 0);
        SGPS_DEM_TRACKED_RESIZE(relPosSphereY, nClumpComponents, "relPosSphereY", 0);
        SGPS_DEM_TRACKED_RESIZE(relPosSphereZ, nClumpComponents, "relPosSphereZ", 0);
    } else {
        SGPS_DEM_TRACKED_RESIZE(radiiSphere, nSpheresGM, "radiiSphere", 0);
        SGPS_DEM_TRACKED_RESIZE(relPosSphereX, nSpheresGM, "relPosSphereX", 0);
        SGPS_DEM_TRACKED_RESIZE(relPosSphereY, nSpheresGM, "relPosSphereY", 0);
        SGPS_DEM_TRACKED_RESIZE(relPosSphereZ, nSpheresGM, "relPosSphereZ", 0);
    }

    // Resize to number of owners
    SGPS_DEM_TRACKED_RESIZE(ownerTypes, nOwnerBodies, "ownerTypes", 0);
    SGPS_DEM_TRACKED_RESIZE(inertiaPropOffsets, nOwnerBodies, "inertiaPropOffsets", 0);
    // If we jitify mass properties, then
    if (solverFlags.useMassJitify) {
        SGPS_DEM_TRACKED_RESIZE(massOwnerBody, nMassProperties, "massOwnerBody", 0);
        SGPS_DEM_TRACKED_RESIZE(mmiXX, nMassProperties, "mmiXX", 0);
        SGPS_DEM_TRACKED_RESIZE(mmiYY, nMassProperties, "mmiYY", 0);
        SGPS_DEM_TRACKED_RESIZE(mmiZZ, nMassProperties, "mmiZZ", 0);
    } else {
        SGPS_DEM_TRACKED_RESIZE(massOwnerBody, nOwnerBodies, "massOwnerBody", 0);
        SGPS_DEM_TRACKED_RESIZE(mmiXX, nOwnerBodies, "mmiXX", 0);
        SGPS_DEM_TRACKED_RESIZE(mmiYY, nOwnerBodies, "mmiYY", 0);
        SGPS_DEM_TRACKED_RESIZE(mmiZZ, nOwnerBodies, "mmiZZ", 0);
    }
    SGPS_DEM_TRACKED_RESIZE(EProxy, nMatTuples, "EProxy", 0);
    SGPS_DEM_TRACKED_RESIZE(nuProxy, nMatTuples, "nuProxy", 0);
    SGPS_DEM_TRACKED_RESIZE(CoRProxy, nMatTuples, "CoRProxy", 0);
    SGPS_DEM_TRACKED_RESIZE(muProxy, nMatTuples, "muProxy", 0);
    SGPS_DEM_TRACKED_RESIZE(CrrProxy, nMatTuples, "CrrProxy", 0);

    // Arrays for contact info
    // The lengths of contact event-based arrays are just estimates. My estimate of total contact pairs is 2n, and I
    // think the max is 6n (although I can't prove it). Note the estimate should be large enough to decrease the number
    // of reallocations in the simulation, but not too large that eats too much memory.
    SGPS_DEM_TRACKED_RESIZE(idGeometryA, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "idGeometryA", 0);
    SGPS_DEM_TRACKED_RESIZE(idGeometryB, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "idGeometryB", 0);
    SGPS_DEM_TRACKED_RESIZE(contactForces, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "contactForces",
                            make_float3(0));
    SGPS_DEM_TRACKED_RESIZE(contactTorque_convToForce, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER,
                            "contactTorque_convToForce", make_float3(0));
    SGPS_DEM_TRACKED_RESIZE(contactHistory, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "contactHistory",
                            make_float3(0));
    SGPS_DEM_TRACKED_RESIZE(contactDuration, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "contactDuration", 0);
    SGPS_DEM_TRACKED_RESIZE(contactType, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "contactType", DEM_NOT_A_CONTACT);
    SGPS_DEM_TRACKED_RESIZE(contactPointGeometryA, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "contactPointGeometryA",
                            make_float3(0));
    SGPS_DEM_TRACKED_RESIZE(contactPointGeometryB, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "contactPointGeometryB",
                            make_float3(0));

    // Transfer buffer arrays
    // The following several arrays will have variable sizes, so here we only used an estimate.
    // It is cudaMalloc-ed memory, not managed, because we want explicit locality control of buffers
    buffer_size = nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER;
    SGPS_DEM_DEVICE_PTR_ALLOC(granData->idGeometryA_buffer, buffer_size);
    SGPS_DEM_DEVICE_PTR_ALLOC(granData->idGeometryB_buffer, buffer_size);
    SGPS_DEM_DEVICE_PTR_ALLOC(granData->contactType_buffer, buffer_size);
    // SGPS_DEM_TRACKED_RESIZE(idGeometryA_buffer, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER, "idGeometryA_buffer",
    // 0); SGPS_DEM_TRACKED_RESIZE(idGeometryB_buffer, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER,
    // "idGeometryB_buffer", 0); SGPS_DEM_TRACKED_RESIZE(contactType_buffer, nOwnerBodies *
    // SGPS_DEM_INIT_CNT_MULTIPLIER, "contactType_buffer", DEM_NOT_A_CONTACT);
    // SGPS_DEM_ADVISE_DEVICE(idGeometryA_buffer, streamInfo.device);
    // SGPS_DEM_ADVISE_DEVICE(idGeometryB_buffer, streamInfo.device);
    // SGPS_DEM_ADVISE_DEVICE(contactType_buffer, streamInfo.device);
    if (!solverFlags.isHistoryless) {
        // SGPS_DEM_TRACKED_RESIZE(contactMapping_buffer, nOwnerBodies * SGPS_DEM_INIT_CNT_MULTIPLIER,
        //                         "contactMapping_buffer", DEM_NULL_MAPPING_PARTNER);
        // SGPS_DEM_ADVISE_DEVICE(contactMapping_buffer, streamInfo.device);
        SGPS_DEM_DEVICE_PTR_ALLOC(granData->contactMapping_buffer, buffer_size);
    }
}

void DEMDynamicThread::registerPolicies(const std::unordered_map<unsigned int, std::string>& template_number_name_map,
                                        const std::vector<float>& clumps_mass_types,
                                        const std::vector<float3>& clumps_moi_types,
                                        const std::vector<float>& ext_obj_mass_types,
                                        const std::vector<float3>& ext_obj_moi_types,
                                        const std::vector<float>& mesh_obj_mass_types,
                                        const std::vector<float3>& mesh_obj_moi_types,
                                        const std::vector<std::shared_ptr<DEMMaterial>>& loaded_materials,
                                        const std::vector<notStupidBool_t>& family_mask_matrix,
                                        const std::set<unsigned int>& no_output_families) {
    // No modification for the arrays in this function. They can only be completely re-constructed.

    // Load in material properties
    for (unsigned int i = 0; i < loaded_materials.size(); i++) {
        std::shared_ptr<DEMMaterial> Mat = loaded_materials.at(i);
        EProxy.at(i) = Mat->E;
        nuProxy.at(i) = Mat->nu;
        CoRProxy.at(i) = Mat->CoR;
        muProxy.at(i) = Mat->mu;
        CrrProxy.at(i) = Mat->Crr;
    }

    // Then load in mass and MOI template info
    size_t k = 0;
    if (solverFlags.useMassJitify) {
        for (unsigned int i = 0; i < clumps_mass_types.size(); i++) {
            massOwnerBody.at(k) = clumps_mass_types.at(i);
            float3 this_moi = clumps_moi_types.at(i);
            mmiXX.at(k) = this_moi.x;
            mmiYY.at(k) = this_moi.y;
            mmiZZ.at(k) = this_moi.z;
            k++;
        }
        for (unsigned int i = 0; i < ext_obj_mass_types.size(); i++) {
            massOwnerBody.at(k) = ext_obj_mass_types.at(i);
            float3 this_moi = ext_obj_moi_types.at(i);
            mmiXX.at(k) = this_moi.x;
            mmiYY.at(k) = this_moi.y;
            mmiZZ.at(k) = this_moi.z;
            k++;
        }
        for (unsigned int i = 0; i < mesh_obj_mass_types.size(); i++) {
            massOwnerBody.at(k) = mesh_obj_mass_types.at(i);
            float3 this_moi = mesh_obj_moi_types.at(i);
            mmiXX.at(k) = this_moi.x;
            mmiYY.at(k) = this_moi.y;
            mmiZZ.at(k) = this_moi.z;
            k++;
        }
    }

    // Store family mask
    for (size_t i = 0; i < family_mask_matrix.size(); i++)
        familyMaskMatrix.at(i) = family_mask_matrix.at(i);

    // Store clump naming map
    templateNumNameMap = template_number_name_map;

    // Take notes of the families that should not be outputted
    {
        std::set<unsigned int>::iterator it;
        unsigned int i = 0;
        familiesNoOutput.resize(no_output_families.size());
        for (it = no_output_families.begin(); it != no_output_families.end(); it++, i++) {
            familiesNoOutput.at(i) = *it;
        }
        std::sort(familiesNoOutput.begin(), familiesNoOutput.end());
        SGPS_DEM_DEBUG_PRINTF("Impl-level families that will not be outputted:");
        SGPS_DEM_DEBUG_EXEC(displayArray<family_t>(familiesNoOutput.data(), familiesNoOutput.size()));
    }
}

void DEMDynamicThread::populateEntityArrays(const std::vector<std::shared_ptr<DEMClumpBatch>>& input_clump_batches,
                                            const std::vector<float3>& input_ext_obj_xyz,
                                            const std::vector<unsigned int>& input_ext_obj_family,
                                            const std::vector<float3>& input_mesh_obj_xyz,
                                            const std::vector<float4>& input_mesh_obj_rot,
                                            const std::vector<unsigned int>& input_mesh_obj_family,
                                            const std::vector<unsigned int>& mesh_facet_owner,
                                            const std::vector<materialsOffset_t>& mesh_facet_materials,
                                            const std::vector<DEMTriangle>& mesh_facets,
                                            const std::vector<std::vector<unsigned int>>& clumps_sp_mat_ids,
                                            const std::vector<float>& clumps_mass_types,
                                            const std::vector<float3>& clumps_moi_types,
                                            const std::vector<std::vector<float>>& clumps_sp_radii_types,
                                            const std::vector<std::vector<float3>>& clumps_sp_location_types,
                                            const std::vector<float>& ext_obj_mass_types,
                                            const std::vector<float3>& ext_obj_moi_types,
                                            const std::vector<float>& mesh_obj_mass_types,
                                            const std::vector<float3>& mesh_obj_moi_types,
                                            size_t nExistOwners,
                                            size_t nExistSpheres) {
    // Load in clump components info (but only if instructed to use jitified clump templates). This step will be
    // repeated even if we are just adding some more clumps to system, not a complete re-initialization.
    size_t k = 0;
    std::vector<unsigned int> prescans_comp;
    if (solverFlags.useClumpJitify) {
        prescans_comp.push_back(0);
        for (const auto& elem : clumps_sp_radii_types) {
            for (const auto& radius : elem) {
                radiiSphere.at(k) = radius;
                k++;
            }
            prescans_comp.push_back(k);
        }
        prescans_comp.pop_back();
        k = 0;

        for (const auto& elem : clumps_sp_location_types) {
            for (const auto& loc : elem) {
                relPosSphereX.at(k) = loc.x;
                relPosSphereY.at(k) = loc.y;
                relPosSphereZ.at(k) = loc.z;
                k++;
            }
        }
    }

    // Left-bottom-front point of the `world'
    float3 LBF;
    LBF.x = simParams->LBFX;
    LBF.y = simParams->LBFY;
    LBF.z = simParams->LBFZ;
    k = 0;
    // Then, load in input clumps
    std::vector<unsigned int> input_clump_types;
    {
        bool pop_family_msg = false;
        std::vector<float3> input_clump_xyz;
        std::vector<float4> input_clump_oriQ;
        std::vector<float3> input_clump_vel;
        std::vector<float3> input_clump_angVel;
        std::vector<unsigned int> input_clump_family;
        // Flatten the input clump batches (because by design we transfer flatten clump info to GPU)
        for (const auto& a_batch : input_clump_batches) {
            // Decode type number and flatten
            std::vector<unsigned int> type_marks(a_batch->GetNumClumps());
            for (size_t i = 0; i < a_batch->GetNumClumps(); i++) {
                type_marks.at(i) = a_batch->types.at(i)->mark;
            }
            input_clump_types.insert(input_clump_types.end(), type_marks.begin(), type_marks.end());
            // Now xyz
            input_clump_xyz.insert(input_clump_xyz.end(), a_batch->xyz.begin(), a_batch->xyz.end());
            // Now vel
            input_clump_vel.insert(input_clump_vel.end(), a_batch->vel.begin(), a_batch->vel.end());
            // Now quaternion
            input_clump_oriQ.insert(input_clump_oriQ.end(), a_batch->oriQ.begin(), a_batch->oriQ.end());
            // Now angular velocity
            input_clump_angVel.insert(input_clump_angVel.end(), a_batch->angVel.begin(), a_batch->angVel.end());
            // For family numbers, we check if the user has explicitly set them. If not, send a warning.
            if ((!a_batch->family_isSpecified) && (!pop_family_msg)) {
                SGPS_DEM_WARNING("Some clumps do not have their family numbers specified, so defaulted to %u",
                                 DEM_DEFAULT_CLUMP_FAMILY_NUM);
                pop_family_msg = true;
            }
            input_clump_family.insert(input_clump_family.end(), a_batch->families.begin(), a_batch->families.end());
            SGPS_DEM_DEBUG_PRINTF("Loaded a batch of %zu clumps.", a_batch->GetNumClumps());
        }

        SGPS_DEM_DEBUG_PRINTF("Total number of transferred clumps this time: %zu", input_clump_types.size());
        SGPS_DEM_DEBUG_PRINTF("Total number of existing owners in simulation: %zu", nExistOwners);
        SGPS_DEM_DEBUG_PRINTF("Total number of owners in simulation after this init call: %zu",
                              simParams->nOwnerBodies);

        for (size_t i = 0; i < input_clump_types.size(); i++) {
            // If got here, this is a clump
            ownerTypes.at(nExistOwners + i) = DEM_OWNER_T_CLUMP;

            auto type_of_this_clump = input_clump_types.at(i);
            inertiaPropOffsets.at(nExistOwners + i) = type_of_this_clump;
            if (!solverFlags.useMassJitify) {
                massOwnerBody.at(nExistOwners + i) = clumps_mass_types.at(type_of_this_clump);
                const float3 this_moi = clumps_moi_types.at(type_of_this_clump);
                mmiXX.at(nExistOwners + i) = this_moi.x;
                mmiYY.at(nExistOwners + i) = this_moi.y;
                mmiZZ.at(nExistOwners + i) = this_moi.z;
            }

            auto this_CoM_coord = input_clump_xyz.at(i) - LBF;

            auto this_clump_no_sp_radii = clumps_sp_radii_types.at(type_of_this_clump);
            auto this_clump_no_sp_relPos = clumps_sp_location_types.at(type_of_this_clump);
            auto this_clump_no_sp_mat_ids = clumps_sp_mat_ids.at(type_of_this_clump);

            for (size_t j = 0; j < this_clump_no_sp_radii.size(); j++) {
                materialTupleOffset.at(nExistSpheres + k) = this_clump_no_sp_mat_ids.at(j);
                ownerClumpBody.at(nExistSpheres + k) = nExistOwners + i;

                // Depending on whether we jitify or flatten
                if (solverFlags.useClumpJitify) {
                    // This component offset, is it too large that can't live in the jitified array?
                    unsigned int this_comp_offset = prescans_comp.at(type_of_this_clump) + j;
                    clumpComponentOffsetExt.at(nExistSpheres + k) = this_comp_offset;
                    if (this_comp_offset < simParams->nJitifiableClumpComponents) {
                        clumpComponentOffset.at(nExistSpheres + k) = this_comp_offset;
                    } else {
                        // If not, an indicator will be put there
                        clumpComponentOffset.at(nExistSpheres + k) = DEM_RESERVED_CLUMP_COMPONENT_OFFSET;
                    }
                } else {
                    radiiSphere.at(nExistSpheres + k) = this_clump_no_sp_radii.at(j);
                    const float3 relPos = this_clump_no_sp_relPos.at(j);
                    relPosSphereX.at(nExistSpheres + k) = relPos.x;
                    relPosSphereY.at(nExistSpheres + k) = relPos.y;
                    relPosSphereZ.at(nExistSpheres + k) = relPos.z;
                }

                k++;
                // std::cout << "Sphere Rel Pos offset: " << this_clump_no_sp_loc_offsets.at(j) << std::endl;
            }

            hostPositionToVoxelID<voxelID_t, subVoxelPos_t, double>(
                voxelID.at(nExistOwners + i), locX.at(nExistOwners + i), locY.at(nExistOwners + i),
                locZ.at(nExistOwners + i), (double)this_CoM_coord.x, (double)this_CoM_coord.y, (double)this_CoM_coord.z,
                simParams->nvXp2, simParams->nvYp2, simParams->voxelSize, simParams->l);

            // Set initial oriQ
            auto oriQ_of_this_clump = input_clump_oriQ.at(i);
            oriQ0.at(nExistOwners + i) = oriQ_of_this_clump.x;
            oriQ1.at(nExistOwners + i) = oriQ_of_this_clump.y;
            oriQ2.at(nExistOwners + i) = oriQ_of_this_clump.z;
            oriQ3.at(nExistOwners + i) = oriQ_of_this_clump.w;

            // Set initial velocity
            auto vel_of_this_clump = input_clump_vel.at(i);
            vX.at(nExistOwners + i) = vel_of_this_clump.x;
            vY.at(nExistOwners + i) = vel_of_this_clump.y;
            vZ.at(nExistOwners + i) = vel_of_this_clump.z;

            // Set initial angular velocity
            auto angVel_of_this_clump = input_clump_angVel.at(i);
            omgBarX.at(nExistOwners + i) = angVel_of_this_clump.x;
            omgBarY.at(nExistOwners + i) = angVel_of_this_clump.y;
            omgBarZ.at(nExistOwners + i) = angVel_of_this_clump.z;

            // Set family code
            family_t this_family_num = input_clump_family.at(i);
            familyID.at(nExistOwners + i) = this_family_num;
        }
    }

    // Load in initial positions and mass properties for the owners of those external objects
    // They go after clump owners
    size_t offset_for_ext_obj = nExistOwners + input_clump_types.size();
    unsigned int offset_for_ext_obj_mass_template = simParams->nDistinctClumpBodyTopologies;
    for (size_t i = 0; i < input_ext_obj_xyz.size(); i++) {
        // If got here, it is an analytical obj
        ownerTypes.at(i + offset_for_ext_obj) = DEM_OWNER_T_ANALYTICAL;
        // Analytical object mass properties are useful in force collection, but not useful in force calculation:
        // analytical component masses are jitified into kernels directly.
        inertiaPropOffsets.at(i + offset_for_ext_obj) = i + offset_for_ext_obj_mass_template;
        if (!solverFlags.useMassJitify) {
            massOwnerBody.at(i + offset_for_ext_obj) = ext_obj_mass_types.at(i);
            const float3 this_moi = ext_obj_moi_types.at(i);
            mmiXX.at(i + offset_for_ext_obj) = this_moi.x;
            mmiYY.at(i + offset_for_ext_obj) = this_moi.y;
            mmiZZ.at(i + offset_for_ext_obj) = this_moi.z;
        }
        auto this_CoM_coord = input_ext_obj_xyz.at(i) - LBF;
        hostPositionToVoxelID<voxelID_t, subVoxelPos_t, double>(
            voxelID.at(i + offset_for_ext_obj), locX.at(i + offset_for_ext_obj), locY.at(i + offset_for_ext_obj),
            locZ.at(i + offset_for_ext_obj), (double)this_CoM_coord.x, (double)this_CoM_coord.y,
            (double)this_CoM_coord.z, simParams->nvXp2, simParams->nvYp2, simParams->voxelSize, simParams->l);
        //// TODO: and initial rot?
        //// TODO: and initial vel?

        family_t this_family_num = input_ext_obj_family.at(i);
        familyID.at(i + offset_for_ext_obj) = this_family_num;
    }

    // Load in initial positions and mass properties for the owners of the meshed objects
    // They go after analytical object owners
    size_t offset_for_mesh_obj = offset_for_ext_obj + input_ext_obj_xyz.size();
    unsigned int offset_for_mesh_obj_mass_template = offset_for_ext_obj_mass_template + input_ext_obj_xyz.size();
    for (size_t i = 0; i < input_mesh_obj_xyz.size(); i++) {
        // If got here, it is a mesh
        ownerTypes.at(i + offset_for_mesh_obj) = DEM_OWNER_T_MESH;

        inertiaPropOffsets.at(i + offset_for_mesh_obj) = i + offset_for_mesh_obj_mass_template;
        if (!solverFlags.useMassJitify) {
            massOwnerBody.at(i + offset_for_mesh_obj) = mesh_obj_mass_types.at(i);
            const float3 this_moi = mesh_obj_moi_types.at(i);
            mmiXX.at(i + offset_for_mesh_obj) = this_moi.x;
            mmiYY.at(i + offset_for_mesh_obj) = this_moi.y;
            mmiZZ.at(i + offset_for_mesh_obj) = this_moi.z;
        }
        auto this_CoM_coord = input_mesh_obj_xyz.at(i) - LBF;
        hostPositionToVoxelID<voxelID_t, subVoxelPos_t, double>(
            voxelID.at(i + offset_for_mesh_obj), locX.at(i + offset_for_mesh_obj), locY.at(i + offset_for_mesh_obj),
            locZ.at(i + offset_for_mesh_obj), (double)this_CoM_coord.x, (double)this_CoM_coord.y,
            (double)this_CoM_coord.z, simParams->nvXp2, simParams->nvYp2, simParams->voxelSize, simParams->l);

        // Set mesh owner's oriQ
        auto oriQ_of_this = input_mesh_obj_rot.at(i);
        oriQ0.at(i + offset_for_mesh_obj) = oriQ_of_this.x;
        oriQ1.at(i + offset_for_mesh_obj) = oriQ_of_this.y;
        oriQ2.at(i + offset_for_mesh_obj) = oriQ_of_this.z;
        oriQ3.at(i + offset_for_mesh_obj) = oriQ_of_this.w;

        //// TODO: and initial vel?

        //// TODO: Per-facet info
        // mesh_facet_owner,
        // mesh_facet_materials,
        // mesh_facets,

        family_t this_family_num = input_mesh_obj_family.at(i);
        familyID.at(i + offset_for_mesh_obj) = this_family_num;
    }
}

void DEMDynamicThread::buildTrackedObjs(const std::vector<std::shared_ptr<DEMClumpBatch>>& input_clump_batches,
                                        const std::vector<float3>& input_ext_obj_xyz,
                                        std::vector<std::shared_ptr<DEMTrackedObj>>& tracked_objs,
                                        size_t nExistOwners) {
    // We take notes on how many clumps each batch has, it will be useful when we assemble the tracker information
    std::vector<size_t> prescans_batch_size;
    prescans_batch_size.push_back(0);
    for (const auto& a_batch : input_clump_batches) {
        prescans_batch_size.push_back(a_batch->GetNumClumps());
    }

    // Provide feedback to the tracked objects, tell them the owner numbers they are looking for
    // Little computation is needed, as long as we know the structure of our owner array: nOwnerClumps go first, then
    // nExtObj, then nTriEntities
    for (auto& tracked_obj : tracked_objs) {
        switch (tracked_obj->type) {
            case (DEM_OWNER_TYPE::CLUMP):
                tracked_obj->ownerID = nExistOwners + prescans_batch_size.at(tracked_obj->load_order);
                break;
            case (DEM_OWNER_TYPE::ANALYTICAL):
                // prescans_batch_size.back() is the total num of loaded clumps this time
                tracked_obj->ownerID = nExistOwners + tracked_obj->load_order + prescans_batch_size.back();
                break;
            default:
                SGPS_DEM_ERROR("A DEM tracked object has an unknown type.");
        }
    }
}

void DEMDynamicThread::initManagedArrays(const std::vector<std::shared_ptr<DEMClumpBatch>>& input_clump_batches,
                                         const std::vector<float3>& input_ext_obj_xyz,
                                         const std::vector<unsigned int>& input_ext_obj_family,
                                         const std::vector<float3>& input_mesh_obj_xyz,
                                         const std::vector<float4>& input_mesh_obj_rot,
                                         const std::vector<unsigned int>& input_mesh_obj_family,
                                         const std::vector<unsigned int>& mesh_facet_owner,
                                         const std::vector<materialsOffset_t>& mesh_facet_materials,
                                         const std::vector<DEMTriangle>& mesh_facets,
                                         const std::unordered_map<unsigned int, std::string>& template_number_name_map,
                                         const std::vector<std::vector<unsigned int>>& clumps_sp_mat_ids,
                                         const std::vector<float>& clumps_mass_types,
                                         const std::vector<float3>& clumps_moi_types,
                                         const std::vector<std::vector<float>>& clumps_sp_radii_types,
                                         const std::vector<std::vector<float3>>& clumps_sp_location_types,
                                         const std::vector<float>& ext_obj_mass_types,
                                         const std::vector<float3>& ext_obj_moi_types,
                                         const std::vector<float>& mesh_obj_mass_types,
                                         const std::vector<float3>& mesh_obj_moi_types,
                                         const std::vector<std::shared_ptr<DEMMaterial>>& loaded_materials,
                                         const std::vector<notStupidBool_t>& family_mask_matrix,
                                         const std::set<unsigned int>& no_output_families,
                                         std::vector<std::shared_ptr<DEMTrackedObj>>& tracked_objs) {
    // Get the info into the managed memory from the host side. Can this process be more efficient? Maybe, but it's
    // initialization anyway.

    registerPolicies(template_number_name_map, clumps_mass_types, clumps_moi_types, ext_obj_mass_types,
                     ext_obj_moi_types, mesh_obj_mass_types, mesh_obj_moi_types, loaded_materials, family_mask_matrix,
                     no_output_families);

    // For initialization, owner array offset is 0
    populateEntityArrays(input_clump_batches, input_ext_obj_xyz, input_ext_obj_family, input_mesh_obj_xyz,
                         input_mesh_obj_rot, input_mesh_obj_family, mesh_facet_owner, mesh_facet_materials, mesh_facets,
                         clumps_sp_mat_ids, clumps_mass_types, clumps_moi_types, clumps_sp_radii_types,
                         clumps_sp_location_types, ext_obj_mass_types, ext_obj_moi_types, mesh_obj_mass_types,
                         mesh_obj_moi_types, 0, 0);

    buildTrackedObjs(input_clump_batches, input_ext_obj_xyz, tracked_objs, 0);
}

void DEMDynamicThread::updateClumpMeshArrays(const std::vector<std::shared_ptr<DEMClumpBatch>>& input_clump_batches,
                                             const std::vector<float3>& input_ext_obj_xyz,
                                             const std::vector<unsigned int>& input_ext_obj_family,
                                             const std::vector<float3>& input_mesh_obj_xyz,
                                             const std::vector<float4>& input_mesh_obj_rot,
                                             const std::vector<unsigned int>& input_mesh_obj_family,
                                             const std::vector<unsigned int>& mesh_facet_owner,
                                             const std::vector<materialsOffset_t>& mesh_facet_materials,
                                             const std::vector<DEMTriangle>& mesh_facets,
                                             const std::vector<std::vector<unsigned int>>& clumps_sp_mat_ids,
                                             const std::vector<float>& clumps_mass_types,
                                             const std::vector<float3>& clumps_moi_types,
                                             const std::vector<std::vector<float>>& clumps_sp_radii_types,
                                             const std::vector<std::vector<float3>>& clumps_sp_location_types,
                                             const std::vector<float>& ext_obj_mass_types,
                                             const std::vector<float3>& ext_obj_moi_types,
                                             const std::vector<float>& mesh_obj_mass_types,
                                             const std::vector<float3>& mesh_obj_moi_types,
                                             const std::vector<std::shared_ptr<DEMMaterial>>& loaded_materials,
                                             const std::vector<notStupidBool_t>& family_mask_matrix,
                                             const std::set<unsigned int>& no_output_families,
                                             std::vector<std::shared_ptr<DEMTrackedObj>>& tracked_objs,
                                             size_t nExistingOwners,
                                             size_t nExistingClumps,
                                             size_t nExistingSpheres,
                                             size_t nExistingTriMesh,
                                             size_t nExistingFacets) {
    // No policy changes here

    // Analytical objects-related arrays should be empty
    populateEntityArrays(input_clump_batches, input_ext_obj_xyz, input_ext_obj_family, input_mesh_obj_xyz,
                         input_mesh_obj_rot, input_mesh_obj_family, mesh_facet_owner, mesh_facet_materials, mesh_facets,
                         clumps_sp_mat_ids, clumps_mass_types, clumps_moi_types, clumps_sp_radii_types,
                         clumps_sp_location_types, ext_obj_mass_types, ext_obj_moi_types, mesh_obj_mass_types,
                         mesh_obj_moi_types, nExistingOwners, nExistingSpheres);

    // No changes to tracked objects
}

void DEMDynamicThread::writeSpheresAsChpf(std::ofstream& ptFile) const {
    chpf::Writer pw;
    // pw.write(ptFile, chpf::Compressor::Type::USE_DEFAULT, mass);
    std::vector<float> posX(simParams->nSpheresGM);
    std::vector<float> posY(simParams->nSpheresGM);
    std::vector<float> posZ(simParams->nSpheresGM);
    std::vector<float> spRadii(simParams->nSpheresGM);
    std::vector<unsigned int> families;
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
        families.resize(simParams->nSpheresGM);
    }
    size_t num_output_spheres = 0;

    for (size_t i = 0; i < simParams->nSpheresGM; i++) {
        auto this_owner = ownerClumpBody.at(i);
        family_t this_family = familyID.at(this_owner);
        // If this (impl-level) family is in the no-output list, skip it
        if (std::binary_search(familiesNoOutput.begin(), familiesNoOutput.end(), this_family)) {
            continue;
        }

        float3 CoM;
        float X, Y, Z;
        voxelID_t voxel = voxelID.at(this_owner);
        subVoxelPos_t subVoxX = locX.at(this_owner);
        subVoxelPos_t subVoxY = locY.at(this_owner);
        subVoxelPos_t subVoxZ = locZ.at(this_owner);
        hostVoxelIDToPosition<float, voxelID_t, subVoxelPos_t>(X, Y, Z, voxel, subVoxX, subVoxY, subVoxZ,
                                                               simParams->nvXp2, simParams->nvYp2, simParams->voxelSize,
                                                               simParams->l);
        CoM.x = X + simParams->LBFX;
        CoM.y = Y + simParams->LBFY;
        CoM.z = Z + simParams->LBFZ;

        size_t compOffset = (solverFlags.useClumpJitify) ? clumpComponentOffsetExt.at(i) : i;
        float this_sp_deviation_x = relPosSphereX.at(compOffset);
        float this_sp_deviation_y = relPosSphereY.at(compOffset);
        float this_sp_deviation_z = relPosSphereZ.at(compOffset);
        float this_sp_rot_0 = oriQ0.at(this_owner);
        float this_sp_rot_1 = oriQ1.at(this_owner);
        float this_sp_rot_2 = oriQ2.at(this_owner);
        float this_sp_rot_3 = oriQ3.at(this_owner);
        hostApplyOriQToVector3<float, float>(this_sp_deviation_x, this_sp_deviation_y, this_sp_deviation_z,
                                             this_sp_rot_0, this_sp_rot_1, this_sp_rot_2, this_sp_rot_3);
        posX.at(num_output_spheres) = CoM.x + this_sp_deviation_x;
        posY.at(num_output_spheres) = CoM.y + this_sp_deviation_y;
        posZ.at(num_output_spheres) = CoM.z + this_sp_deviation_z;
        // std::cout << "Sphere Pos: " << posX.at(i) << ", " << posY.at(i) << ", " << posZ.at(i) << std::endl;

        spRadii.at(num_output_spheres) = radiiSphere.at(compOffset);

        // Family number
        if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
            families.at(num_output_spheres) = this_family;
        }

        num_output_spheres++;
    }
    // Write basics
    posX.resize(num_output_spheres);
    posY.resize(num_output_spheres);
    posZ.resize(num_output_spheres);
    spRadii.resize(num_output_spheres);
    // TODO: Set {} to the list of column names
    pw.write(ptFile, chpf::Compressor::Type::USE_DEFAULT,
             {DEM_OUTPUT_FILE_X_COL_NAME, DEM_OUTPUT_FILE_Y_COL_NAME, DEM_OUTPUT_FILE_Z_COL_NAME, "r"}, posX, posY,
             posZ, spRadii);
    // Write family numbers
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
        families.resize(num_output_spheres);
        // TODO: How to do that?
        // pw.write(ptFile, chpf::Compressor::Type::USE_DEFAULT, {}, families);
    }
}

void DEMDynamicThread::writeSpheresAsCsv(std::ofstream& ptFile) const {
    std::ostringstream outstrstream;

    outstrstream << DEM_OUTPUT_FILE_X_COL_NAME + "," + DEM_OUTPUT_FILE_Y_COL_NAME + "," + DEM_OUTPUT_FILE_Z_COL_NAME +
                        ",r";
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ABSV) {
        outstrstream << ",absv";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::VEL) {
        outstrstream << ",v_x,v_y,v_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ANG_VEL) {
        outstrstream << ",w_x,w_y,w_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ACC) {
        outstrstream << ",a_x,a_y,a_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ANG_ACC) {
        outstrstream << ",alpha_x,alpha_y,alpha_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
        outstrstream << ",family";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::MAT) {
        outstrstream << ",material";
    }
    outstrstream << "\n";

    for (size_t i = 0; i < simParams->nSpheresGM; i++) {
        auto this_owner = ownerClumpBody.at(i);
        family_t this_family = familyID.at(this_owner);
        // If this (impl-level) family is in the no-output list, skip it
        if (std::binary_search(familiesNoOutput.begin(), familiesNoOutput.end(), this_family)) {
            continue;
        }

        float3 CoM;
        float3 pos;
        float radius;
        float X, Y, Z;
        voxelID_t voxel = voxelID.at(this_owner);
        subVoxelPos_t subVoxX = locX.at(this_owner);
        subVoxelPos_t subVoxY = locY.at(this_owner);
        subVoxelPos_t subVoxZ = locZ.at(this_owner);
        hostVoxelIDToPosition<float, voxelID_t, subVoxelPos_t>(X, Y, Z, voxel, subVoxX, subVoxY, subVoxZ,
                                                               simParams->nvXp2, simParams->nvYp2, simParams->voxelSize,
                                                               simParams->l);
        CoM.x = X + simParams->LBFX;
        CoM.y = Y + simParams->LBFY;
        CoM.z = Z + simParams->LBFZ;

        size_t compOffset = (solverFlags.useClumpJitify) ? clumpComponentOffsetExt.at(i) : i;
        float3 this_sp_deviation;
        this_sp_deviation.x = relPosSphereX.at(compOffset);
        this_sp_deviation.y = relPosSphereY.at(compOffset);
        this_sp_deviation.z = relPosSphereZ.at(compOffset);
        float this_sp_rot_0 = oriQ0.at(this_owner);
        float this_sp_rot_1 = oriQ1.at(this_owner);
        float this_sp_rot_2 = oriQ2.at(this_owner);
        float this_sp_rot_3 = oriQ3.at(this_owner);
        hostApplyOriQToVector3<float, float>(this_sp_deviation.x, this_sp_deviation.y, this_sp_deviation.z,
                                             this_sp_rot_0, this_sp_rot_1, this_sp_rot_2, this_sp_rot_3);
        pos = CoM + this_sp_deviation;
        outstrstream << pos.x << "," << pos.y << "," << pos.z;

        radius = radiiSphere.at(compOffset);
        outstrstream << "," << radius;

        // Only linear velocity
        if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ABSV) {
            float3 absv;
            absv.x = vX.at(this_owner);
            absv.y = vY.at(this_owner);
            absv.z = vZ.at(this_owner);
            outstrstream << "," << length(absv);
        }

        // Family number needs to be user number
        if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
            outstrstream << "," << +(this_family);
        }

        outstrstream << "\n";
    }

    ptFile << outstrstream.str();
}

void DEMDynamicThread::writeClumpsAsChpf(std::ofstream& ptFile) const {}

void DEMDynamicThread::writeClumpsAsCsv(std::ofstream& ptFile) const {
    std::ostringstream outstrstream;

    // xyz and quaternion are always there
    outstrstream << DEM_OUTPUT_FILE_X_COL_NAME + "," + DEM_OUTPUT_FILE_Y_COL_NAME + "," + DEM_OUTPUT_FILE_Z_COL_NAME +
                        ",Q0,Q1,Q2,Q3," + DEM_OUTPUT_FILE_CLUMP_TYPE_NAME;
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ABSV) {
        outstrstream << ",absv";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::VEL) {
        outstrstream << ",v_x,v_y,v_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ANG_VEL) {
        outstrstream << ",w_x,w_y,w_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ACC) {
        outstrstream << ",a_x,a_y,a_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ANG_ACC) {
        outstrstream << ",alpha_x,alpha_y,alpha_z";
    }
    if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
        outstrstream << ",family";
    }
    outstrstream << "\n";

    for (size_t i = 0; i < simParams->nOwnerBodies; i++) {
        // i is this owner's number. And if it is not a clump, we can move on.
        if (ownerTypes.at(i) != DEM_OWNER_T_CLUMP)
            continue;

        family_t this_family = familyID.at(i);
        // If this (impl-level) family is in the no-output list, skip it
        if (std::binary_search(familiesNoOutput.begin(), familiesNoOutput.end(), this_family)) {
            continue;
        }

        float3 CoM;
        float X, Y, Z;
        voxelID_t voxel = voxelID.at(i);
        subVoxelPos_t subVoxX = locX.at(i);
        subVoxelPos_t subVoxY = locY.at(i);
        subVoxelPos_t subVoxZ = locZ.at(i);
        hostVoxelIDToPosition<float, voxelID_t, subVoxelPos_t>(X, Y, Z, voxel, subVoxX, subVoxY, subVoxZ,
                                                               simParams->nvXp2, simParams->nvYp2, simParams->voxelSize,
                                                               simParams->l);
        CoM.x = X + simParams->LBFX;
        CoM.y = Y + simParams->LBFY;
        CoM.z = Z + simParams->LBFZ;
        // Output position
        outstrstream << CoM.x << "," << CoM.y << "," << CoM.z;

        // Then quaternions
        outstrstream << "," << oriQ0.at(i) << "," << oriQ1.at(i) << "," << oriQ2.at(i) << "," << oriQ3.at(i);

        // Then type of clump
        unsigned int clump_mark = inertiaPropOffsets.at(i);
        outstrstream << "," << templateNumNameMap.at(clump_mark);

        // Only linear velocity
        if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::ABSV) {
            float3 absv;
            absv.x = vX.at(i);
            absv.y = vY.at(i);
            absv.z = vZ.at(i);
            outstrstream << "," << length(absv);
        }

        // Family number needs to be user number
        if (solverFlags.outputFlags & DEM_OUTPUT_CONTENT::FAMILY) {
            outstrstream << "," << +(this_family);
        }

        outstrstream << "\n";
    }

    ptFile << outstrstream.str();
}

inline void DEMDynamicThread::contactEventArraysResize(size_t nContactPairs) {
    SGPS_DEM_TRACKED_RESIZE_NOPRINT(idGeometryA, nContactPairs);
    SGPS_DEM_TRACKED_RESIZE_NOPRINT(idGeometryB, nContactPairs);
    SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactType, nContactPairs);
    SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactForces, nContactPairs);
    SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactTorque_convToForce, nContactPairs);

    SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactPointGeometryA, nContactPairs);
    SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactPointGeometryB, nContactPairs);

    // Re-pack pointers in case the arrays got reallocated
    granData->idGeometryA = idGeometryA.data();
    granData->idGeometryB = idGeometryB.data();
    granData->contactType = contactType.data();
    granData->contactForces = contactForces.data();
    granData->contactTorque_convToForce = contactTorque_convToForce.data();
    granData->contactPointGeometryA = contactPointGeometryA.data();
    granData->contactPointGeometryB = contactPointGeometryB.data();
}

inline void DEMDynamicThread::unpackMyBuffer() {
    // Make a note on the contact number of the previous time step
    *stateOfSolver_resources.pNumPrevContacts = *stateOfSolver_resources.pNumContacts;

    GPU_CALL(cudaMemcpy(stateOfSolver_resources.pNumContacts, &(granData->nContactPairs_buffer), sizeof(size_t),
                        cudaMemcpyDeviceToDevice));

    // Need to resize those contact event-based arrays before usage
    if (*stateOfSolver_resources.pNumContacts > idGeometryA.size() ||
        *stateOfSolver_resources.pNumContacts > buffer_size) {
        contactEventArraysResize(*stateOfSolver_resources.pNumContacts);
    }

    GPU_CALL(cudaMemcpy(granData->idGeometryA, granData->idGeometryA_buffer,
                        *stateOfSolver_resources.pNumContacts * sizeof(bodyID_t), cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->idGeometryB, granData->idGeometryB_buffer,
                        *stateOfSolver_resources.pNumContacts * sizeof(bodyID_t), cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->contactType, granData->contactType_buffer,
                        *stateOfSolver_resources.pNumContacts * sizeof(contact_t), cudaMemcpyDeviceToDevice));
    if (!solverFlags.isHistoryless) {
        // Note we don't have to use dedicated memory space for unpacking contactMapping_buffer contents, because we
        // only use it once per kT update, at the time of unpacking. So let us just use a temp vector to store it. Note
        // we cannot use vector 0 since it may hold critical flattened owner ID info.
        size_t mapping_bytes = (*stateOfSolver_resources.pNumContacts) * sizeof(contactPairs_t);
        granData->contactMapping = (contactPairs_t*)stateOfSolver_resources.allocateTempVector(1, mapping_bytes);
        GPU_CALL(cudaMemcpy(granData->contactMapping, granData->contactMapping_buffer, mapping_bytes,
                            cudaMemcpyDeviceToDevice));
    }
}

inline void DEMDynamicThread::sendToTheirBuffer() {
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_voxelID, granData->voxelID,
                        simParams->nOwnerBodies * sizeof(voxelID_t), cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_locX, granData->locX, simParams->nOwnerBodies * sizeof(subVoxelPos_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_locY, granData->locY, simParams->nOwnerBodies * sizeof(subVoxelPos_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_locZ, granData->locZ, simParams->nOwnerBodies * sizeof(subVoxelPos_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ0, granData->oriQ0, simParams->nOwnerBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ1, granData->oriQ1, simParams->nOwnerBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ2, granData->oriQ2, simParams->nOwnerBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_oriQ3, granData->oriQ3, simParams->nOwnerBodies * sizeof(oriQ_t),
                        cudaMemcpyDeviceToDevice));

    // Family number is a typical changable quantity on-the-fly. If this flag is on, dT is responsible for sending this
    // info to kT.
    if (solverFlags.canFamilyChange) {
        GPU_CALL(cudaMemcpy(granData->pKTOwnedBuffer_familyID, granData->familyID,
                            simParams->nOwnerBodies * sizeof(family_t), cudaMemcpyDeviceToDevice));
    }
}

inline void DEMDynamicThread::migrateContactHistory() {
    // Use this newHistory and newDuration to store temporarily the rearranged contact history.  Note we cannot use
    // vector 0 or 1 since they may be in use.
    size_t history_arr_bytes = (*stateOfSolver_resources.pNumContacts) * sizeof(float3);
    size_t duration_arr_bytes = (*stateOfSolver_resources.pNumContacts) * sizeof(float);
    size_t sentry_bytes = (*stateOfSolver_resources.pNumPrevContacts) * sizeof(notStupidBool_t);
    float3* newHistory = (float3*)stateOfSolver_resources.allocateTempVector(2, history_arr_bytes);
    float* newDuration = (float*)stateOfSolver_resources.allocateTempVector(3, duration_arr_bytes);
    // This is used for checking if there are contact history got lost in the transition by surprise. But no need to
    // check if the user did not ask for it.
    notStupidBool_t* contactSentry = (notStupidBool_t*)stateOfSolver_resources.allocateTempVector(4, sentry_bytes);

    // A sentry array is here to see if there exist a contact that dT thinks it's alive but kT doesn't map it to the new
    // history array
    size_t blocks_needed_for_rearrange;
    if (*stateOfSolver_resources.pNumPrevContacts > 0 && verbosity >= DEM_VERBOSITY::STEP_METRIC) {
        // GPU_CALL(cudaMemset(contactSentry, 0, sentry_bytes));
        blocks_needed_for_rearrange = (*stateOfSolver_resources.pNumPrevContacts + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) /
                                      SGPS_DEM_MAX_THREADS_PER_BLOCK;
        if (blocks_needed_for_rearrange > 0) {
            prep_force_kernels->kernel("markAliveContacts")
                .instantiate()
                .configure(dim3(blocks_needed_for_rearrange), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0,
                           streamInfo.stream)
                .launch(granData->contactDuration, contactSentry, *stateOfSolver_resources.pNumPrevContacts);
            GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
        }
    }

    // Rearrange contact histories based on kT instruction
    blocks_needed_for_rearrange =
        (*stateOfSolver_resources.pNumContacts + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) / SGPS_DEM_MAX_THREADS_PER_BLOCK;
    if (blocks_needed_for_rearrange > 0) {
        prep_force_kernels->kernel("rearrangeContactHistory")
            .instantiate()
            .configure(dim3(blocks_needed_for_rearrange), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0, streamInfo.stream)
            .launch(granData->contactMapping, granData->contactHistory, granData->contactDuration, newHistory,
                    newDuration, contactSentry, *stateOfSolver_resources.pNumContacts);
        GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
    }

    // Take a look, does the sentry indicate that there is an `alive' contact got lost?
    if (*stateOfSolver_resources.pNumPrevContacts > 0 && verbosity >= DEM_VERBOSITY::STEP_METRIC) {
        notStupidBool_t* lostContact =
            (notStupidBool_t*)stateOfSolver_resources.allocateTempVector(5, sizeof(notStupidBool_t));
        boolMaxReduce(contactSentry, lostContact, *stateOfSolver_resources.pNumPrevContacts, streamInfo.stream,
                      stateOfSolver_resources);
        if (*lostContact && solverFlags.isAsync) {
            SGPS_DEM_STEP_METRIC(
                "At least one contact is active at time %.9g on dT, but it is not detected on kT, therefore being "
                "removed unexpectedly!",
                timeElapsed);
            SGPS_DEM_DEBUG_PRINTF("New contact A:");
            SGPS_DEM_DEBUG_EXEC(displayArray<bodyID_t>(granData->idGeometryA, *stateOfSolver_resources.pNumContacts));
            SGPS_DEM_DEBUG_PRINTF("New contact B:");
            SGPS_DEM_DEBUG_EXEC(displayArray<bodyID_t>(granData->idGeometryB, *stateOfSolver_resources.pNumContacts));
            SGPS_DEM_DEBUG_PRINTF("Old contact history:");
            SGPS_DEM_DEBUG_EXEC(displayFloat3(granData->contactHistory, *stateOfSolver_resources.pNumPrevContacts));
            SGPS_DEM_DEBUG_PRINTF("History mapping:");
            SGPS_DEM_DEBUG_EXEC(
                displayArray<contactPairs_t>(granData->contactMapping, *stateOfSolver_resources.pNumContacts));
            SGPS_DEM_DEBUG_PRINTF("Sentry:");
            SGPS_DEM_DEBUG_EXEC(
                displayArray<notStupidBool_t>(contactSentry, *stateOfSolver_resources.pNumPrevContacts));
        }
    }

    // Copy new history back to history array (after resizing the `main' history array)
    if (*stateOfSolver_resources.pNumContacts > contactHistory.size()) {
        SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactHistory, *stateOfSolver_resources.pNumContacts);
        SGPS_DEM_TRACKED_RESIZE_NOPRINT(contactDuration, *stateOfSolver_resources.pNumContacts);
        granData->contactHistory = contactHistory.data();
        granData->contactDuration = contactDuration.data();
    }
    GPU_CALL(cudaMemcpy(granData->contactHistory, newHistory, (*stateOfSolver_resources.pNumContacts) * sizeof(float3),
                        cudaMemcpyDeviceToDevice));
    GPU_CALL(cudaMemcpy(granData->contactDuration, newDuration, (*stateOfSolver_resources.pNumContacts) * sizeof(float),
                        cudaMemcpyDeviceToDevice));
}

inline void DEMDynamicThread::calculateForces() {
    // reset force (acceleration) arrays for this time step and apply gravity
    size_t threads_needed_for_prep = simParams->nOwnerBodies > *stateOfSolver_resources.pNumContacts
                                         ? simParams->nOwnerBodies
                                         : *stateOfSolver_resources.pNumContacts;
    size_t blocks_needed_for_prep =
        (threads_needed_for_prep + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) / SGPS_DEM_MAX_THREADS_PER_BLOCK;

    prep_force_kernels->kernel("prepareForceArrays")
        .instantiate()
        .configure(dim3(blocks_needed_for_prep), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0, streamInfo.stream)
        .launch(simParams, granData, *stateOfSolver_resources.pNumContacts);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));

    // TODO: is there a better way??? Like memset?
    // GPU_CALL(cudaMemset(granData->contactForces, zeros, *stateOfSolver_resources.pNumContacts * sizeof(float3)));
    // GPU_CALL(cudaMemset(granData->alphaX, 0, simParams->nOwnerBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->alphaY, 0, simParams->nOwnerBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->alphaZ, 0, simParams->nOwnerBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->aX,
    //                     (double)simParams->h * (double)simParams->h * (double)simParams->Gx / (double)simParams->l,
    //                     simParams->nOwnerBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->aY,
    //                     (double)simParams->h * (double)simParams->h * (double)simParams->Gy / (double)simParams->l,
    //                     simParams->nOwnerBodies * sizeof(float)));
    // GPU_CALL(cudaMemset(granData->aZ,
    //                     (double)simParams->h * (double)simParams->h * (double)simParams->Gz / (double)simParams->l,
    //                     simParams->nOwnerBodies * sizeof(float)));

    size_t blocks_needed_for_contacts =
        (*stateOfSolver_resources.pNumContacts + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) / SGPS_DEM_NUM_BODIES_PER_BLOCK;
    // If no contact then we don't have to calculate forces. Note there might still be forces, coming from prescription
    // or other sources.
    if (blocks_needed_for_contacts > 0) {
        timers.GetTimer("Calculate contact forces").start();
        // a custom kernel to compute forces
        cal_force_kernels->kernel("calculateContactForces")
            .instantiate()
            .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, streamInfo.stream)
            .launch(simParams, granData, *stateOfSolver_resources.pNumContacts);
        GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
        // displayFloat3(granData->contactForces, *stateOfSolver_resources.pNumContacts);
        // std::cout << "===========================" << std::endl;
        timers.GetTimer("Calculate contact forces").stop();

        timers.GetTimer("Collect contact forces").start();
        // Reflect those body-wise forces on their owner clumps
        // hostCollectForces(granData->inertiaPropOffsets, granData->idGeometryA, granData->idGeometryB,
        //                   granData->contactForces, granData->aX, granData->aY, granData->aZ,
        //                   granData->ownerClumpBody, granData->massOwnerBody, simParams->h,
        //                   *stateOfSolver_resources.pNumContacts,simParams->l);
        collectContactForces(collect_force_kernels, granData, *stateOfSolver_resources.pNumContacts,
                             simParams->nOwnerBodies, contactPairArr_isFresh, streamInfo.stream,
                             stateOfSolver_resources, timers);
        // displayArray<float>(granData->aX, simParams->nOwnerBodies);
        // displayFloat3(granData->contactForces, *stateOfSolver_resources.pNumContacts);
        // std::cout << *stateOfSolver_resources.pNumContacts << std::endl;

        // Calculate the torque on owner clumps from those body-wise forces
        // hostCollectTorques(granData->inertiaPropOffsets, granData->idGeometryA, granData->idGeometryB,
        //                    granData->contactForces, granData->contactPointGeometryA, granData->contactPointGeometryB,
        //                    granData->alphaX, granData->alphaY, granData->alphaZ, granData->ownerClumpBody,
        //                    granData->mmiXX, granData->mmiYY, granData->mmiZZ, simParams->h,
        //                    *stateOfSolver_resources.pNumContacts, simParams->l);
        // displayArray<float>(granData->oriQ0, simParams->nOwnerBodies);
        // displayArray<float>(granData->oriQ1, simParams->nOwnerBodies);
        // displayArray<float>(granData->oriQ2, simParams->nOwnerBodies);
        // displayArray<float>(granData->oriQ3, simParams->nOwnerBodies);
        // displayArray<float>(granData->alphaX, simParams->nOwnerBodies);
        // displayArray<float>(granData->alphaY, simParams->nOwnerBodies);
        // displayArray<float>(granData->alphaZ, simParams->nOwnerBodies);
        timers.GetTimer("Collect contact forces").stop();
    }
}

inline void DEMDynamicThread::integrateClumpMotions() {
    size_t blocks_needed_for_clumps =
        (simParams->nOwnerBodies + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) / SGPS_DEM_NUM_BODIES_PER_BLOCK;
    integrator_kernels->kernel("integrateClumps")
        .instantiate()
        .configure(dim3(blocks_needed_for_clumps), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, streamInfo.stream)
        .launch(granData, simParams->nOwnerBodies, simParams->h, timeElapsed);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
}

inline void DEMDynamicThread::routineChecks() {
    if (solverFlags.canFamilyChange) {
        size_t blocks_needed_for_clumps =
            (simParams->nOwnerBodies + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) / SGPS_DEM_NUM_BODIES_PER_BLOCK;
        mod_kernels->kernel("applyFamilyChanges")
            .instantiate()
            .configure(dim3(blocks_needed_for_clumps), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, streamInfo.stream)
            .launch(granData, simParams->nOwnerBodies, simParams->h, timeElapsed);
        GPU_CALL(cudaStreamSynchronize(streamInfo.stream));
    }
}

void DEMDynamicThread::workerThread() {
    // Set the gpu for this thread
    GPU_CALL(cudaSetDevice(streamInfo.device));
    GPU_CALL(cudaStreamCreate(&streamInfo.stream));

    while (!pSchedSupport->dynamicShouldJoin) {
        {
            std::unique_lock<std::mutex> lock(pSchedSupport->dynamicStartLock);
            while (!pSchedSupport->dynamicStarted) {
                pSchedSupport->cv_DynamicStartLock.wait(lock);
            }
            // Ensure that we wait for start signal on next iteration
            pSchedSupport->dynamicStarted = false;
            // The following is executed when kT and dT are being destroyed
            if (pSchedSupport->dynamicShouldJoin) {
                break;
            }
        }

        // There is only one situation where dT needs to wait for kT to provide one initial CD result...
        // This is the `new-boot' case, where stampLastUpdateOfDynamic == -1; in any other situations, dT does not have
        // `drift-into-future-too-much' problem here, b/c if it has the problem then it would have been addressed at the
        // end of last DoDynamics call, the final `ShouldWait' check.
        if (pSchedSupport->stampLastUpdateOfDynamic < 0) {
            // In this `new-boot' case, we send kT a work order, b/c dT needs results from CD to proceed. After this one
            // instance, kT and dT may work in an async fashion.
            {
                std::lock_guard<std::mutex> lock(pSchedSupport->kinematicOwnedBuffer_AccessCoordination);
                sendToTheirBuffer();
            }
            pSchedSupport->kinematicOwned_Cons2ProdBuffer_isFresh = true;
            contactPairArr_isFresh = true;
            pSchedSupport->schedulingStats.nKinematicUpdates++;
            // Signal the kinematic that it has data for a new work order.
            pSchedSupport->cv_KinematicCanProceed.notify_all();
            // Then dT will wait for kT to finish one initial run
            {
                std::unique_lock<std::mutex> lock(pSchedSupport->dynamicCanProceed);
                while (!pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh) {
                    // loop to avoid spurious wakeups
                    pSchedSupport->cv_DynamicCanProceed.wait(lock);
                }
            }
        }

        for (double cycle = 0.0; cycle < cycleDuration; cycle += simParams->h) {
            // If the produce is fresh, use it
            if (pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh) {
                timers.GetTimer("Unpack updates from kT").start();
                {
                    // Acquire lock and use the content of the dynamic-owned transfer buffer
                    std::lock_guard<std::mutex> lock(pSchedSupport->dynamicOwnedBuffer_AccessCoordination);
                    unpackMyBuffer();
                    // Leave myself a mental note that I just obtained new produce from kT
                    contactPairArr_isFresh = true;
                    // pSchedSupport->schedulingStats.nDynamicReceives++;
                }
                // dT got the produce, now mark its buffer to be no longer fresh
                pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh = false;
                pSchedSupport->stampLastUpdateOfDynamic = (pSchedSupport->currentStampOfDynamic).load();

                // If this is a history-based run, then when contacts are received, we need to migrate the contact
                // history info, to match the structure of the new contact array
                if (!solverFlags.isHistoryless) {
                    migrateContactHistory();
                }
                timers.GetTimer("Unpack updates from kT").stop();
            }

            // If using variable ts size, only when a step is accepted can we move on
            bool step_accepted = false;
            do {
                calculateForces();

                routineChecks();

                timers.GetTimer("Integration").start();
                integrateClumpMotions();
                timers.GetTimer("Integration").stop();

                step_accepted = true;
            } while ((!solverFlags.isStepConst) || (!step_accepted));

            // CalculateForces is done, set contactPairArr_isFresh to false
            // This will be set to true next time it receives an update from kT
            contactPairArr_isFresh = false;

            /*
            if (cycle == (cycleDuration - 1))
                pSchedSupport->dynamicDone = true;
            */

            // Dynamic wrapped up one cycle, record this fact into schedule support
            pSchedSupport->currentStampOfDynamic++;

            // If the kinematic is idle, give it the opportunity to get busy again
            if (!pSchedSupport->kinematicOwned_Cons2ProdBuffer_isFresh) {
                timers.GetTimer("Send to kT buffer").start();
                // Acquire lock and refresh the work order for the kinematic
                {
                    std::lock_guard<std::mutex> lock(pSchedSupport->kinematicOwnedBuffer_AccessCoordination);
                    sendToTheirBuffer();
                }
                pSchedSupport->kinematicOwned_Cons2ProdBuffer_isFresh = true;
                pSchedSupport->schedulingStats.nKinematicUpdates++;
                timers.GetTimer("Send to kT buffer").stop();
                // Signal the kinematic that it has data for a new work order
                pSchedSupport->cv_KinematicCanProceed.notify_all();
            }

            // Check if we need to wait; i.e., if dynamic drifted too much into future, then we must wait a bit before
            // the next cycle begins
            if (pSchedSupport->dynamicShouldWait()) {
                timers.GetTimer("Wait for kT update").start();
                // Wait for a signal from kT to indicate that kT has caught up
                pSchedSupport->schedulingStats.nTimesDynamicHeldBack++;
                std::unique_lock<std::mutex> lock(pSchedSupport->dynamicCanProceed);
                while (!pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh) {
                    // Loop to avoid spurious wakeups
                    pSchedSupport->cv_DynamicCanProceed.wait(lock);
                }
                timers.GetTimer("Wait for kT update").stop();
            }
            // NOTE: This ShouldWait check is at the end of a dT cycle, not at the beginning, because dT could start a
            // cycle and immediately be too much into future. If this is at the beginning, dT will stale while kT could
            // also be chilling waiting for an update; but since it's at the end of a cycle, kT should got busy already
            // due to that kinematicOwned_Cons2ProdBuffer_isFresh got set to true just before.

            // TODO: make changes for variable time step size cases
            timeElapsed += simParams->h;

            nTotalSteps++;
        }

        // When getting here, dT has finished one user call (although perhaps not at the end of the user script)
        pPagerToMain->userCallDone = true;
        pPagerToMain->cv_mainCanProceed.notify_all();
    }
}

void DEMDynamicThread::getTiming(std::vector<std::string>& names, std::vector<double>& vals) {
    names = timer_names;
    for (const auto& name : timer_names) {
        vals.push_back(timers.GetTimer(name).GetTimeSeconds());
    }
}

void DEMDynamicThread::startThread() {
    std::lock_guard<std::mutex> lock(pSchedSupport->dynamicStartLock);
    pSchedSupport->dynamicStarted = true;
    pSchedSupport->cv_DynamicStartLock.notify_one();
}

void DEMDynamicThread::resetUserCallStat() {
    // Reset last kT-side data receiving cycle time stamp.
    pSchedSupport->stampLastUpdateOfDynamic = -1;
    pSchedSupport->currentStampOfDynamic = 0;
    // Reset dT stats variables, making ready for next user call
    pSchedSupport->dynamicDone = false;
    pSchedSupport->dynamicOwned_Prod2ConsBuffer_isFresh = false;
    contactPairArr_isFresh = true;
}

size_t DEMDynamicThread::estimateMemUsage() const {
    return m_approx_bytes_used;
}

void DEMDynamicThread::jitifyKernels(const std::unordered_map<std::string, std::string>& Subs) {
    // First one is force array preparation kernels
    {
        prep_force_kernels = std::make_shared<jitify::Program>(
            std::move(JitHelper::buildProgram("DEMPrepForceKernels", JitHelper::KERNEL_DIR / "DEMPrepForceKernels.cu",
                                              Subs, {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
    }
    // Then force calculation kernels
    // Depending on historyless-ness, we may want jitify different versions of kernels
    {
        if (solverFlags.isHistoryless) {
            cal_force_kernels = std::make_shared<jitify::Program>(std::move(JitHelper::buildProgram(
                "DEMHistorylessForceKernels", JitHelper::KERNEL_DIR / "DEMHistorylessForceKernels.cu", Subs,
                {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
        } else {
            cal_force_kernels = std::make_shared<jitify::Program>(std::move(JitHelper::buildProgram(
                "DEMHistoryBasedForceKernels", JitHelper::KERNEL_DIR / "DEMHistoryBasedForceKernels.cu", Subs,
                {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
        }
    }
    // Then force accumulation kernels
    {
        collect_force_kernels = std::make_shared<jitify::Program>(std::move(
            JitHelper::buildProgram("DEMCollectForceKernels", JitHelper::KERNEL_DIR / "DEMCollectForceKernels.cu", Subs,
                                    {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
    }
    // Then integration kernels
    {
        integrator_kernels = std::make_shared<jitify::Program>(std::move(
            JitHelper::buildProgram("DEMIntegrationKernels", JitHelper::KERNEL_DIR / "DEMIntegrationKernels.cu", Subs,
                                    {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
    }
    // Then kernels that are... wildcards, which make on-the-fly changes to solver data
    if (solverFlags.canFamilyChange) {
        mod_kernels = std::make_shared<jitify::Program>(
            std::move(JitHelper::buildProgram("DEMModeratorKernels", JitHelper::KERNEL_DIR / "DEMModeratorKernels.cu",
                                              Subs, {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
    }
    // Then quarrying kernels
    // {
    //     quarry_stats_kernels = std::make_shared<jitify::Program>(
    //         std::move(JitHelper::buildProgram("DEMQueryKernels", JitHelper::KERNEL_DIR / "DEMQueryKernels.cu", Subs,
    //                                           {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
    // }
    // Then misc kernels
    {
        misc_kernels = std::make_shared<jitify::Program>(
            std::move(JitHelper::buildProgram("DEMMiscKernels", JitHelper::KERNEL_DIR / "DEMMiscKernels.cu", Subs,
                                              {"-I" + (JitHelper::KERNEL_DIR / "..").string()})));
    }
}

float* DEMDynamicThread::inspectCall(const std::shared_ptr<jitify::Program>& inspection_kernel,
                                     const std::string& kernel_name,
                                     size_t n,
                                     DEM_CUB_REDUCE_FLAVOR reduce_flavor,
                                     bool all_domain) {
    // We can use temp vectors as we please
    size_t quarryTempSize = n * sizeof(float);
    float* resArr = (float*)stateOfSolver_resources.allocateTempVector(1, quarryTempSize);
    size_t regionTempSize = n * sizeof(notStupidBool_t);
    notStupidBool_t* boolArr = (notStupidBool_t*)stateOfSolver_resources.allocateTempVector(2, regionTempSize);
    GPU_CALL(cudaMemset(boolArr, 0, regionTempSize));

    size_t returnSize = sizeof(float);
    float* res = (float*)stateOfSolver_resources.allocateTempVector(3, returnSize);
    size_t blocks_needed = (n + SGPS_DEM_MAX_THREADS_PER_BLOCK - 1) / SGPS_DEM_MAX_THREADS_PER_BLOCK;
    inspection_kernel->kernel(kernel_name)
        .instantiate()
        .configure(dim3(blocks_needed), dim3(SGPS_DEM_MAX_THREADS_PER_BLOCK), 0, streamInfo.stream)
        .launch(granData, resArr, boolArr, n);
    GPU_CALL(cudaStreamSynchronize(streamInfo.stream));

    if (all_domain) {
        switch (reduce_flavor) {
            case (DEM_CUB_REDUCE_FLAVOR::MAX):
                floatMaxReduce(resArr, res, n, streamInfo.stream, stateOfSolver_resources);
                break;
            case (DEM_CUB_REDUCE_FLAVOR::SUM):
                floatSumReduce(resArr, res, n, streamInfo.stream, stateOfSolver_resources);
                break;
            case (DEM_CUB_REDUCE_FLAVOR::NONE):
                return resArr;
        }
        // If this inspection is comfined in a region, then boolArr and resArr need to be sorted and reduce by key
    } else {
        //// TODO: Implement it
    }

    return res;
}

float3 DEMDynamicThread::getOwnerAngVel(bodyID_t ownerID) const {
    float3 angVel;
    angVel.x = omgBarX.at(ownerID);
    angVel.y = omgBarY.at(ownerID);
    angVel.z = omgBarZ.at(ownerID);
    return angVel;
}

float4 DEMDynamicThread::getOwnerOriQ(bodyID_t ownerID) const {
    float4 oriQ;
    oriQ.x = oriQ0.at(ownerID);
    oriQ.y = oriQ1.at(ownerID);
    oriQ.z = oriQ2.at(ownerID);
    oriQ.w = oriQ3.at(ownerID);
    return oriQ;
}

float3 DEMDynamicThread::getOwnerVel(bodyID_t ownerID) const {
    float3 vel;
    vel.x = vX.at(ownerID);
    vel.y = vY.at(ownerID);
    vel.z = vZ.at(ownerID);
    return vel;
}

float3 DEMDynamicThread::getOwnerPos(bodyID_t ownerID) const {
    float3 pos;
    double X, Y, Z;
    voxelID_t voxel = voxelID.at(ownerID);
    subVoxelPos_t subVoxX = locX.at(ownerID);
    subVoxelPos_t subVoxY = locY.at(ownerID);
    subVoxelPos_t subVoxZ = locZ.at(ownerID);
    hostVoxelIDToPosition<double, voxelID_t, subVoxelPos_t>(X, Y, Z, voxel, subVoxX, subVoxY, subVoxZ, simParams->nvXp2,
                                                            simParams->nvYp2, simParams->voxelSize, simParams->l);
    pos.x = X + simParams->LBFX;
    pos.y = Y + simParams->LBFY;
    pos.z = Z + simParams->LBFZ;
    return pos;
}

void DEMDynamicThread::setOwnerAngVel(bodyID_t ownerID, float3 angVel) {
    omgBarX.at(ownerID) = angVel.x;
    omgBarY.at(ownerID) = angVel.y;
    omgBarZ.at(ownerID) = angVel.z;
}

void DEMDynamicThread::setOwnerPos(bodyID_t ownerID, float3 pos) {
    // Convert to relative pos wrt LBF point first
    double X, Y, Z;
    X = pos.x - simParams->LBFX;
    Y = pos.y - simParams->LBFY;
    Z = pos.z - simParams->LBFZ;
    hostPositionToVoxelID<voxelID_t, subVoxelPos_t, double>(voxelID.at(ownerID), locX.at(ownerID), locY.at(ownerID),
                                                            locZ.at(ownerID), X, Y, Z, simParams->nvXp2,
                                                            simParams->nvYp2, simParams->voxelSize, simParams->l);
}

void DEMDynamicThread::setOwnerOriQ(bodyID_t ownerID, float4 oriQ) {
    oriQ0.at(ownerID) = oriQ.x;
    oriQ1.at(ownerID) = oriQ.y;
    oriQ2.at(ownerID) = oriQ.z;
    oriQ3.at(ownerID) = oriQ.w;
}

void DEMDynamicThread::setOwnerVel(bodyID_t ownerID, float3 vel) {
    vX.at(ownerID) = vel.x;
    vY.at(ownerID) = vel.y;
    vZ.at(ownerID) = vel.z;
}

}  // namespace sgps
