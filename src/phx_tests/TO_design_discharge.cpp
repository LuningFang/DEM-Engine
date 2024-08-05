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

    int particle_family = 0;
    int plate_family = 2;
    int recylcled_family = 1;

    if (argc != 3){
        std::cout << "Usage: ./TO_design_setup <Test ID> <flow rate num, 0 slow, 1 fast>" << std::endl;
        return 1;
    }

    std::string TO_design_mesh;
    std::string TEST_NAME;
    std::string out_dir;
    std::string orifice_filename;


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
    std::string input_particle_positions = TEST_NAME + "/settling/settled.csv";


    int flow_rate = std::stoi(argv[2]);
    if (flow_rate == 0){
        out_dir = TEST_NAME + "/discharge_slow/";
        orifice_filename = "clumps/TO_bottom_plate_20_percent.csv";
    }
    else if (flow_rate == 1){
        out_dir = TEST_NAME + "/discharge_fast/";
        orifice_filename = "clumps/TO_bottom_plate_50_percent.csv";

    }
    else {
        std::cout << "Usage: ./TO_design_setup <Test ID, 0 for uniform and 1 for allConstr> <flow rate num, 0 slow, 1 fast>" << std::endl;
        return 1;
    }
    create_directories(out_dir);


    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL | OUTPUT_CONTENT::FAMILY);
    DEMSim.SetErrorOutVelocity(1e7);

    float fric_coef = 0.6;
    float carbo_density = 3.6;
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    std::vector<double> radius_array = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};

    float step_size = 2.5e-6;
    float time_end = 5.0;
    unsigned int fps = 100;

    // extract test name from TO_design_mesh
    // need to put move this somewhere else ... 
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

    std::cout << "added TO unit cell " << std::endl; 

    // Add sim domain
    DEMSim.InstructBoxDomainDimension({ 0.,          bxDim}, 
                                      { 0.,          byDim},
                                      {-bzDim     , bzDim/  2. + 5});
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_wall);

    // Load particles from a csv file
    LoadParticlesFromFile(DEMSim, input_particle_positions, mat_type_carbo, radius_array, carbo_density, particle_family);

    ///////////////// Add orifice /////////////////////
    AddOrificeParticles(DEMSim, orifice_filename, mat_type_wall, plate_family);

    // set up periodic boundary, where particles discharged added to the top
    DEMSim.SetFamilyPrescribedPosition(recylcled_family, "none", "none", "Z-22");
    DEMSim.SetFamilyPrescribedLinVel(recylcled_family, "0", "0", "none");

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, 981));
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.Initialize();


    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // write the mesh file to outdir 
    DEMSim.WriteMeshFile(out_dir + "/TO_design.vtk");

    std::pair<float, float> plate_x_range = std::make_pair(0, bxDim);
    std::pair<float, float> plate_y_range = std::make_pair(0, byDim);
    std::pair<float, float> plate_z_range = std::make_pair(18., 20.);

    for (double t = 0; t < (double)time_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            DEMSim.ShowThreadCollaborationStats();
            
            std::chrono::high_resolution_clock::time_point start_step = std::chrono::high_resolution_clock::now();

            // change family of particles that are discharged
            int num_changed = DEMSim.ChangeClumpFamily(recylcled_family, plate_x_range, plate_y_range, plate_z_range);
            std::cout << "number of family changed: " << num_changed << std::endl;
            DEMSim.DoDynamicsThenSync(0.);

            std::cout << "Frame: " << currframe << std::endl;
            std::cout << "Time: " << t << std::endl;
            float recycled_family_mass = DEMSim.GetFamilyMass(recylcled_family);
            std::cout << "Recycled family mass: " << recycled_family_mass << " g, mass flow rate: " <<  recycled_family_mass / (1./fps) << "g / sec"   << std::endl;

            char filename[200];
            sprintf(filename, "%s/particles_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            std::cout << "write file: " << filename << std::endl;
            currframe++;
            std::chrono::high_resolution_clock::time_point end_step = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end_step - start_step);
            std::cout << "Time for writing and family changing: " << time_sec.count() << " seconds" << std::endl;

        }

        DEMSim.DoDynamics(step_size);
        if (curr_step % out_steps == 0) {
            DEMSim.ChangeFamily(recylcled_family, particle_family);
            DEMSim.DoDynamicsThenSync(0.);

        }
    }


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count())  << " seconds (wall time) to finish " << time_end << " seconds simulation"
              << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    return 0;
}