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

    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;

    if (argc != 2){
        std::cout << "Usage: ./phx_setup <Test ID>" << std::endl;
        return 1;
    }

    int test_id = std::stoi(argv[1]);

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY);
    DEMSim.SetErrorOutVelocity(1e7);
    DEMSim.SetNoForceRecord();

    float fric_coef = 0.6;
    float carbo_density = 3.6;
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    std::vector<double> radius_array = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};

    float step_size = 2e-6;
    float time_end = 2.0;
    unsigned int fps = 10;

    std::string TEST_NAME = "DEBUG_" + std::to_string(test_id);

    auto mat_type_carbo = DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", fric_coef}, {"Crr", 0.0}});
    auto mat_type_wall = DEMSim.LoadMaterial({{"E", 2e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", fric_coef}, {"Crr", 0.0}});

    std::string pin_pos_file = "sim_data/pin_pos_run" + std::to_string(test_id) + ".csv";
    std::vector<float3> pin_centers;
    double pin_hdim;
    // mesh pins
    if (test_id != 0){
        std::string pin_mesh_file_name = "mesh/pins/teardrop_run" + std::to_string(test_id) + ".obj";
        pin_hdim = 0.45;
        // Add mesh pins
        pin_centers = AddMeshPins(DEMSim, pin_pos_file, pin_mesh_file_name, mat_type_wall);
    }
    else {
        // cylindrical pins
        pin_hdim = 0.2;
        pin_centers = AddCylindricalPins(DEMSim, pin_pos_file, pin_hdim, mat_type_wall);
    }

    // Add sim domain
    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.}, 
                                      {-byDim / 2., byDim},
                                      {-bzDim / 2., bzDim/  2.});
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_wall);


    // Add particles
    AddParticles(DEMSim, radius_array, carbo_density, mat_type_carbo, make_float3(bxDim, byDim, bzDim), pin_centers, pin_hdim);


    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, -981, 0));
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/" + TEST_NAME + "/settling/";
    create_directories(out_dir);

    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

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

int AddParticles(DEMSolver& DEMSim, 
                  std::vector<double> radius_array, 
                  float density, 
                  std::shared_ptr<DEMMaterial> mat_type, 
                  float3 box_dim,
                  std::vector<float3> pin_centers,
                  float pin_hdim) {

    // Create clump templates for the filling materials
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    // Then randomly create some clump templates for filling the drum    
    for (int i = 0; i < radius_array.size(); i++) {
        double volume = 4./3. * PI * std::pow(radius_array[i],3);
        // Load a sphere type
        clump_types.push_back(DEMSim.LoadSphereType(volume * density, radius_array[i], mat_type));
    }

    // Then sample some particles inside the drum
    std::vector<std::shared_ptr<DEMClumpTemplate>> particle_template_type;
    std::vector<float3> particle_xyz;

    box_dim.x = box_dim.x/2.;
    box_dim.y = box_dim.y/8.;
    box_dim.z = box_dim.z/2.;
    particle_xyz =  PopulateParticlePositions(2.02 * radius_array[0], box_dim, pin_centers, pin_hdim);
    std::cout << "number of particles: " << particle_xyz.size() << std::endl;

    unsigned int num_clumps = particle_xyz.size();
    // TODO: this should be changed based on the size ratio 
    for (unsigned int i = 0; i < num_clumps; i++) {
        particle_template_type.push_back(clump_types.at(i % clump_types.size()));
    }
    DEMSim.AddClumps(particle_template_type, particle_xyz);
    return num_clumps;
}

