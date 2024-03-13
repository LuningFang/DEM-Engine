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


int main() {

    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY | OUTPUT_CONTENT::VEL);
    // DEMSim.SetVerbosity(STEP_METRIC);

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    DEMSim.SetNoForceRecord();

    auto mat_type_sand = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.5}, {"Crr", 0.01}});
    auto mat_type_plate = DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.5}, {"Crr", 0.01}});


    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.}, 
                                      {-36.       , byDim / 2.},
                                      {-bzDim / 2., bzDim/  2.});


    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_plate);


    // what does top open do? 
    // DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Since two types of materials have the same mu, this following call does not change the default mu for their
    // interaction, it's still 0.5.
    DEMSim.SetMaterialPropertyPair("mu", mat_type_sand, mat_type_plate, 0.5);

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
        pin->AddCylinder(pin_center, CylAxis, tube_radius, mat_type_plate, 1); // outward normal
        std::cout << "added cylinder at " << pin_center.x << ", " << pin_center.y << ", " << pin_center.z << std::endl;
    }

    // load the particles


    auto part1_clump_xyz = DEMSim.ReadClumpXyzFromCsv("./DemoOutput_phx_pins/settled.csv");
    auto part1_clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./DemoOutput_phx_pins/settled.csv");

    std::vector<float3> in_xyz;
    std::vector<float4> in_quat;
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    unsigned int t_num = 0;
    unsigned int sphere_types = 3;
    for (int i = 0; i < sphere_types; i++) {
        char t_name[20];
        sprintf(t_name, "%04d", t_num);

        auto this_type_xyz = part1_clump_xyz[std::string(t_name)];
        auto this_type_quat = part1_clump_quaternion[std::string(t_name)];

        size_t n_clump_this_type = this_type_xyz.size();
        // Prepare clump type identification vector for loading into the system (don't forget type 0 in
        // ground_particle_templates is the template for rover wheel)
        std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                 clump_types.at(t_num));

        // Add them to the big long vector
        in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
        in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
        in_types.insert(in_types.end(), this_type.begin(), this_type.end());

        // Our template names are 0000, 0001 etc.
        t_num++;
    }
    // Finally, load them into the system
    DEMClumpBatch base_batch(in_xyz.size());
    base_batch.SetTypes(in_types);
    base_batch.SetPos(in_xyz);
    base_batch.SetOriQ(in_quat);
    base_batch.SetFamily(0);
    DEMSim.AddClumps(base_batch);


    // Load bottom plate particles

    // Read plate particles from data file in data/sim_data/plate_pos.csv, skip the first header row, and assemble the value into a list of float3
    std::vector<float3> plate_particles;

    std::ifstream plate_pos_file(GetDEMEDataFile("clumps/bottom_plate.csv"));
    // file has 4 columns, x, y, z, and r, and the first row is the header
    // we only need the first 3 columns for the position
    std::string line;
    std::getline(plate_pos_file, line); // skip the first line
    while (std::getline(plate_pos_file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        float3 particle_pos = make_float3(std::stof(tokens[0]), std::stof(tokens[1]), std::stof(tokens[2]));
        plate_particles.push_back(particle_pos);
    }

    auto plate_template =
        DEMSim.LoadClumpType(2, make_float3(0.5, 0.5, 0.5),
                             std::vector<float>(plate_particles.size(), 0.025), plate_particles, mat_type_plate);
    std::cout << plate_particles.size() << " spheres make up the rotating drum" << std::endl;

    // Add drum
    auto plate = DEMSim.AddClumps(plate_template, make_float3(0));
    unsigned int plate_family = 100;
    plate->SetFamilies(plate_family);

    DEMSim.SetFamilyFixed(plate_family);
    // Disable contacts within drum components
    DEMSim.DisableContactBetweenFamilies(plate_family, plate_family);

////////////////////////////////////////////////////




    // Keep tab of the max velocity in simulation
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");
    float max_v;

    float step_size = 1e-6;
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
    out_dir += "/DemoOutput_phx_discharge";
    create_directory(out_dir);

    float time_end = 1.5;
    unsigned int fps = 2000;
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
