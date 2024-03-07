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


// Read pin positions from a file
std::vector<float3> ReadPinPositions(std::string filename) {
    std::vector<float3> pin_positions;

    std::ifstream pin_pos_file(GetDEMEDataFile(filename));
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
        pin_positions.push_back(pin_center);
    }
    return pin_positions;
}


std::vector<float3> PopulateParticlePositions(float spacing, 
                                              float3 dim, 
                                              std::vector<float3> pin_pos, 
                                              double pin_radius) {

    std::vector<float3> material_points;

    // Generate initial clumps for piling
    PDSampler sampler(spacing); // spacing should be 2.01 * particle_radius 

    float box_X = dim.x;
    float box_Y = dim.y;
    float box_Z = dim.z;

    float fill_y_bottom = -box_Y / 2.0f + 2 * spacing;
    float fill_y_top = box_Y / 2.0f - spacing;
    float3 sampler_center = make_float3(0.0f, fill_y_bottom, 0.0f);
    float3 sampler_hdim = make_float3(box_X / 2.0f - spacing, 1e-6, box_Z / 2.0f - spacing);

    // iterate over the pin position provided in pin_pos, find the unique values of pin_pos.y and keep them in an array of float
    std::vector<float> unique_y_pos;
    for (auto pin : pin_pos) {
        double y = pin.y;
        // check to see if the y value is already in the unique_y_pos array
        if (std::find(unique_y_pos.begin(), unique_y_pos.end(), y) == unique_y_pos.end()) {
            unique_y_pos.push_back(y);
        }
    }

    std::sort(unique_y_pos.begin(), unique_y_pos.end());
    int tube_pos_itr = 0;  // iterator over unique_y_pos, start layer number from bottom

    while (sampler_center.y <= fill_y_top) {

        if (tube_pos_itr >= unique_y_pos.size()){
            break;
        }

        // check if sample_center.y() belongs to their range of the unique y position of the tube
        bool in_tube = false;

        double curr_tube_lower_bound = unique_y_pos[tube_pos_itr] - pin_radius;
        double curr_tube_upper_bound = unique_y_pos[tube_pos_itr] + pin_radius;

        if (sampler_center.y < curr_tube_lower_bound - 0.5 * spacing) {
            // we are in no tube zone, create sample points however we want
            auto points = sampler.SampleBox(sampler_center, sampler_hdim);
            material_points.insert(material_points.end(), points.begin(), points.end());
            std::cout << "insert points at layer " << sampler_center.y << std::endl;
            sampler_center.y += spacing;
        } else {
            // oh no we hit the tube zone
            // we need to create sample points that are not inside the tube
            sampler_center.y = curr_tube_upper_bound + 0.5 * spacing;
            tube_pos_itr++;
        }
    }

    // now we are at the top of the tube, keep adding layers until we reach the top of the box
    while (sampler_center.y <= fill_y_top) {
        auto points = sampler.SampleBox(sampler_center, sampler_hdim);
        material_points.insert(material_points.end(), points.begin(), points.end());
        std::cout << "insert more points at layer " << sampler_center.y << std::endl;
        sampler_center.y += spacing;
    }

    return material_points;
}


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
    double radius_array[3] = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};

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
    std::vector<float3> pin_centers = ReadPinPositions("sim_data/pin_pos.csv");
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

    input_xyz =  PopulateParticlePositions(2.02 * radius_array[0], make_float3(bxDim, byDim, bzDim), pin_centers, tube_radius);

    std::cout << "number of particles: " << input_xyz.size() << std::endl;

    unsigned int num_clumps = input_xyz.size();
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

    char cp_filename[200];
    sprintf(cp_filename, "%s/GRC_3e5.csv", out_dir.c_str());
    DEMSim.WriteClumpFile(std::string(cp_filename));

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    char cnt_filename[200];
    sprintf(cnt_filename, "%s/Contact_pairs_3e5.csv", out_dir.c_str());
    DEMSim.WriteContactFile(std::string(cnt_filename));


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
