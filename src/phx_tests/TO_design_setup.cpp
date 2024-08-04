//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// TODO: at the end of the script, print out stats like bulk density, 
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>
#include <heat_exchanger_util.h>

#include <cstdio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

// given particle radius array, density, material type, bounding box, generate random particle positions, and add particles to the simulation given random radius size 
int AddParticles(DEMSolver& DEMSim, 
                  std::vector<double> radius_array, 
                  float density, 
                  std::shared_ptr<DEMMaterial> mat_type, 
                  float3 box_dim,
                  std::vector<float3> pin_centers,
                  float pin_hdim);

int main(int argc, char** argv) {

    double bxDim = 0.5;
    double byDim = 5.0; // change this to 5 for the actual test 
    double bzDim = 34;

    if (argc != 2){
        std::cout << "Usage: ./TO_design_setup <Test ID>" << std::endl;
        return 1;
    }

    std::string TO_design_mesh;
    std::string TEST_NAME;

    if (std::string(argv[1]) == "0") {
        TO_design_mesh = "mesh/TO/uniform.obj";
        TEST_NAME = "TO_uniform";
    } else if (std::string(argv[1]) == "1") {
        TO_design_mesh = "mesh/TO/allConstr.obj";   
        TEST_NAME = "TO_allConstr";
    } else {
        std::cout << "Usage: ./TO_design_setup <Test ID, 0 for uniform and 1 for allConstr>" << std::endl;
        return 1;
    }

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY);
    DEMSim.SetErrorOutVelocity(1e7);
    DEMSim.SetNoForceRecord();

    float fric_coef = 0.6;
    float carbo_density = 3.6;
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    std::vector<double> radius_array = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};

    float step_size = 2.5e-6;
    float time_end = 5.0;
    unsigned int fps = 100;

    // extract test name from TO_design_mesh
    auto mat_type_carbo = DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", fric_coef}, {"Crr", 0.0}});
    auto mat_type_wall = DEMSim.LoadMaterial({{"E", 2e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", fric_coef}, {"Crr", 0.0}});

    auto TO_mesh = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile(TO_design_mesh), mat_type_wall);
    std::cout << TO_mesh->GetNumTriangles() << " faces and " << TO_mesh->GetNumNodes() << " vertices" << std::endl;

    TO_mesh = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile(TO_design_mesh), mat_type_wall);
    TO_mesh->Mirror(make_float3(0, 1, 0), make_float3(0, 1, 0));

    TO_mesh = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile(TO_design_mesh), mat_type_wall);
    TO_mesh->Move(make_float3(0,2,0), make_float4(0, 0, 0, 1));

    TO_mesh = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile(TO_design_mesh), mat_type_wall);
    TO_mesh->Mirror(make_float3(0, 1, 0), make_float3(0, 1, 0));
    TO_mesh->Move(make_float3(0,2,0), make_float4(0, 0, 0, 1));

    TO_mesh = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile(TO_design_mesh), mat_type_wall);
    TO_mesh->Move(make_float3(0,4,0), make_float4(0, 0, 0, 1));



    // pin->Move(pin_center, rot);
    // pin->SetFamily(family_ID);
    int family_ID = 10;
    std::cout << "added TO unit cell " << std::endl; 
    DEMSim.SetFamilyFixed(family_ID);

    // Add sim domain
    DEMSim.InstructBoxDomainDimension({ 0.,          bxDim}, 
                                      { 0.,          byDim},
                                      {-bzDim     , bzDim/  2.});
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_wall);


    // Create clump templates for the filling materials
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Then randomly create some clump templates for filling the drum    
    for (int i = 0; i < radius_array.size(); i++) {
        double volume = 4./3. * PI * std::pow(radius_array[i],3);
        // Load a sphere type
        clump_types.push_back(DEMSim.LoadSphereType(volume * carbo_density, radius_array[i], mat_type_carbo));
    }

    // Then sample some particles inside the drum
    std::vector<std::shared_ptr<DEMClumpTemplate>> particle_template_type;
    std::vector<float3> particle_xyz;

    // sampler spacing
    PDSampler smapler(2.02 * radius_array[0]);
    float3 sampler_center = make_float3(bxDim/2., byDim/2., -bzDim/2.);
    float3 sampler_hdim = make_float3(bxDim/2. * 0.95, byDim/2. * 0.95, bzDim/2. * 0.9);
    particle_xyz = smapler.SampleBox(sampler_center, sampler_hdim);

    unsigned int num_clumps = particle_xyz.size();
    // TODO: this should be changed based on the size ratio 
    for (unsigned int i = 0; i < num_clumps; i++) {
        particle_template_type.push_back(clump_types.at(i % clump_types.size()));
    }
    DEMSim.AddClumps(particle_template_type, particle_xyz);    


    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, 981));
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/" + TEST_NAME + "/settling/";
    create_directories(out_dir);

    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();


    // write the mesh file to outdir 
    DEMSim.WriteMeshFile(out_dir.string() + "TO_design.vtk");


    for (double t = 0; t < (double)time_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[100];
            sprintf(filename, "%s/particles_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            currframe++;
        }
        DEMSim.DoDynamics(step_size);
    }

    char cp_filename[200];
    sprintf(cp_filename, "%s/settled.csv", out_dir.c_str());
    DEMSim.WriteClumpFile(std::string(cp_filename));


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count())  << " seconds (wall time) to finish " << time_end << " seconds simulation"
              << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    return 0;
}