//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// This demo features an analytical boundary-represented fast rotating container
// with particles of various shapes pulled into it. Different types of particles
// are marked with different family numbers (identification numbers) for easier
// visualizations.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;


int main() {

    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY);
    // DEMSim.SetVerbosity(STEP_METRIC);

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    DEMSim.SetNoForceRecord();

    auto mat_type_sand = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.5}, {"Crr", 0.01}});
    auto mat_type_drum = DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.5}, {"Crr", 0.01}});


    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.}, 
                                      {-byDim / 2., byDim / 2.},
                                      {-bzDim / 2., bzDim/  2.});


    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_drum);


    // what does top open do? 
    // DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Since two types of materials have the same mu, this following call does not change the default mu for their
    // interaction, it's still 0.5.
    DEMSim.SetMaterialPropertyPair("mu", mat_type_sand, mat_type_drum, 0.5);

    // Create some random clump templates for the filling materials
    // An array to store these generated clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    double sand_density = 2.6e3;
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    double radius_array[3] = {0.212 * scaling, 0.15 * scaling, 0.106 * scaling};

    // ratio, 20%, 50%, 30%
    double mass;

    // Then randomly create some clump templates for filling the drum    
    for (int i = 0; i < 3; i++) {
        double volume = 4./3. * PI * std::pow(radius_array[i],3);
        mass = volume * sand_density;
        // Load a sphere type
        clump_types.push_back(DEMSim.LoadSphereType(mass, radius_array[i], mat_type_sand));
    }

    // Add the cylinder pins
    float3 CylAxis = make_float3(0, 0, 1);
    float tube_radius = 0.2;
    float CylHeight = bzDim*1.2;
    float safe_delta = 0.03;

    // read pin center from data file in data/sim_data/pin_pos.csv, skip the first header row, and assemble the value into a list of float3
    std::vector<float3> pin_centers;
    std::ifstream pin_pos_file(GetDEMEDataFile("sim_data/pin_pos.csv"));
    std::string line;
    std::getline(pin_pos_file, line); // skip the first line
    while (std::getline(pin_pos_file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        float3 pin_center = make_float3(std::stof(tokens[0]), std::stof(tokens[1]), std::stof(tokens[2]));
        pin_centers.push_back(pin_center);
    }

    std::cout << "number of pins: " << pin_centers.size() << std::endl;

    for (auto pin_center : pin_centers) {
        auto pin = DEMSim.AddExternalObject();
        pin->AddCylinder(pin_center, CylAxis, tube_radius, mat_type_drum, 1); // outward normal
        std::cout << "added cylinder at " << pin_center.x << ", " << pin_center.y << ", " << pin_center.z << std::endl;
    }
    // Then sample some particles inside the drum
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_template_type;
    std::vector<float3> input_xyz;
    std::vector<unsigned int> family_code;

    float3 sample_center = make_float3(0, byDim/4., 0);

    double smaple_half_dim_y = byDim/4. * 0.95; // <-- coefficient set to 0.2 for testing, 0.8 for actual running
    float3 sample_half_dim = make_float3(bxDim/2. * 0.95, smaple_half_dim_y, bzDim/2. * 0.95);

    std::cout << "sample center: " << sample_center.x << ", " << sample_center.y << ", " << sample_center.z << std::endl;
    std::cout << "sample h_dim: " << sample_half_dim.x << ", " << sample_half_dim.y << ", " << sample_half_dim.z << std::endl;

    auto input_material_xyz =
        DEMBoxGridSampler(sample_center, sample_half_dim,
                          radius_array[0] * 2.1,  radius_array[0] * 2.1, radius_array[0] * 2.1);
    input_xyz.insert(input_xyz.end(), input_material_xyz.begin(), input_material_xyz.end());
    unsigned int num_clumps = input_material_xyz.size();
    // Casually select from generated clump types
    for (unsigned int i = 0; i < num_clumps; i++) {
        input_template_type.push_back(clump_types.at(i % clump_types.size()));
        // Every clump type that has a unique mass, gets a unique family number
        family_code.push_back((i % clump_types.size()) / 2);
    }
    auto particles = DEMSim.AddClumps(input_template_type, input_xyz);
    particles->SetFamilies(family_code);

    // Keep tab of the max velocity in simulation
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");
    float max_v;

    float step_size = 5e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, -981, 0));
    DEMSim.SetExpandSafetyType("auto");
    // If there is a velocity that an analytical object (i.e. the drum) has that you'd like the solver to take into
    // account in consideration of adding contact margins, you have to specify it here, since the solver's automatic max
    // velocity derivation algorithm currently cannot take analytical object's angular velocity-induced velocity into
    // account.
    DEMSim.SetExpandSafetyAdder(6.0);
    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/DemoOutput_phx_pins";
    create_directory(out_dir);

    float time_end = 3.0;
    unsigned int fps = 20;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (double t = 0; t < (double)time_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[100];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            currframe++;
            max_v = max_v_finder->GetValue();
            std::cout << "Max velocity of any point in simulation is " << max_v << std::endl;

        }

        DEMSim.DoDynamics(step_size);
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count()) / time_end * 10.0 << " seconds (wall time) to finish 10 seconds' simulation"
              << std::endl;
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    std::cout << "DEMdemo_Centrifuge exiting..." << std::endl;
    return 0;
}
