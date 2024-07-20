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
#include <heat_exchanger_util.h>

#include <cstdio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;


int main() {

    int recycled_family = 10;


    // read the structure from file 
    std::vector<float3> tpms_positions = ReadTPMSPositions("clumps/tpms_unit_cell.csv");

    std::cout << "number of tpms particles in one unit cell: " << tpms_positions.size() << std::endl;

    // now find the bounding box of the structure
    float3 min = make_float3(1e6, 1e6, 1e6);
    float3 max = make_float3(-1e6, -1e6, -1e6);
    for (auto pos : tpms_positions) {
        if (pos.x < min.x) min.x = pos.x;
        if (pos.y < min.y) min.y = pos.y;
        if (pos.z < min.z) min.z = pos.z;
        if (pos.x > max.x) max.x = pos.x;
        if (pos.y > max.y) max.y = pos.y;
        if (pos.z > max.z) max.z = pos.z;
    }

    double tpms_diameter = 0.08;
    double tpms_radius = tpms_diameter / 2.0;
    double particle_density = 2.5;

    // unit cell dimension
    double unit_x = max.x - min.x + tpms_diameter;
    double unit_y = max.y - min.y + tpms_diameter;
    double unit_z = max.z - min.z + tpms_diameter;


    int3 unit_cell_num = make_int3(3,1,1); 

    // Domain should be dependent on the TPMS size and the number of unit cells in x,y and z directions
    double bxDim = unit_x * unit_cell_num.x;
    double byDim = unit_y * unit_cell_num.y;
    double bzDim = unit_z * unit_cell_num.z;

    std::cout << "unit cell dimension: " << unit_x << " " << unit_y << " " << unit_z << std::endl;
    std::cout << "tpms structure size: " << bxDim << " " << byDim << " " << bzDim << std::endl;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY | OUTPUT_CONTENT::VEL);
    // DEMSim.SetVerbosity(STEP_METRIC);

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    DEMSim.SetNoForceRecord();

    auto mat_type_sand = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.16}, {"Crr", 0.01}});
    auto mat_type_tpms = DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.16}, {"Crr", 0.01}});


    // direction of the flow is y, so make byDim larger than tpms structure
    DEMSim.InstructBoxDomainDimension(bxDim, 5 * byDim, bzDim);
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_tpms);


    DEMSim.SetMaterialPropertyPair("mu", mat_type_sand, mat_type_tpms, 0.5);

    // Create some random clump templates for the filling materials
    // An array to store these generated clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    double sand_density = 2.6e3;
    double scaling = 0.1;  // for testing, actual particle scale is 0.1

    double radius_array[4] = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling, 0.125 * scaling};

    // ratio, 20%, 50%, 30%
    double mass;

    // Then randomly create some clump templates for filling the drum    
    for (int i = 0; i < 4; i++) {
        double volume = 4./3. * PI * std::pow(radius_array[i],3);
        mass = volume * sand_density;
        // Load a sphere type
        clump_types.push_back(DEMSim.LoadSphereType(mass, radius_array[i], mat_type_sand));
    }



    // Sample particles above the structure
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_template_type;
    std::vector<float3> input_xyz;
    std::vector<unsigned int> family_code;


    // define sample particle center and dimensions
    double spacing = 2.01 * radius_array[0];
    PDSampler sampler(spacing);
    float fill_y_bottom = byDim/2. + spacing;
    float fill_y_top = 2.5 * byDim - spacing;
    float3 sample_center = make_float3(0, fill_y_bottom, 0);
    float3 sample_hdim = make_float3(bxDim / 2. - spacing, 1e-6, bzDim / 2. - spacing);

    while (sample_center.y < fill_y_top) {
        auto points = sampler.SampleBox(sample_center, sample_hdim);
        input_xyz.insert(input_xyz.end(), points.begin(), points.end());
        sample_center.y += spacing;
    }

    unsigned int num_clumps = input_xyz.size();
    // Casually select from generated clump types
    for (unsigned int i = 0; i < num_clumps; i++) {
        input_template_type.push_back(clump_types.at(i % clump_types.size()));
        // Every clump type that has a unique mass, gets a unique family number
        family_code.push_back((i % clump_types.size()));
    }
    auto particles = DEMSim.AddClumps(input_template_type, input_xyz);
    particles->SetFamilies(family_code);


    //////////////////////////////////////////
    //// Load tpms points first, then add ////
    //////////////////////////////////////////
    std::vector<float3> tpms_particles;

    std::ifstream tpms_file(GetDEMEDataFile("clumps/tpms_tiled_density_50.csv"));
    // file has 4 columns, x, y, z, and r, and the first row is the header
    // we only need the first 3 columns for the position
    std::string line;
    std::getline(tpms_file, line); // skip the first line
    while (std::getline(tpms_file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }
        float3 pos = make_float3(std::stof(tokens[0]), std::stof(tokens[1]), std::stof(tokens[2]));
        tpms_particles.push_back(pos);
        // tpms_particles.push_back(make_float3(pos.x + unit_x, pos.y, pos.z));
        // tpms_particles.push_back(make_float3(pos.x - unit_x, pos.y, pos.z));
    }

    auto tpms_template =
        DEMSim.LoadClumpType(2, make_float3(0.5, 0.5, 0.5),
                             std::vector<float>(tpms_particles.size(), tpms_radius), tpms_particles, mat_type_tpms);
    std::cout << tpms_particles.size() << " spheres make up the rotating drum" << std::endl;

    // write tpms particles to a file
    std::ofstream tpms_out("tpms_particles_output.csv");
    tpms_out << "x,y,z,r" << std::endl;
    for (auto pos : tpms_particles) {
        tpms_out << pos.x << "," << pos.y << "," << pos.z << "," << tpms_radius << std::endl;
    }


    // Add tpms
    auto tpms = DEMSim.AddClumps(tpms_template, make_float3(0));
    unsigned int tpms_family = 100;
    tpms->SetFamilies(tpms_family);

    

    DEMSim.SetFamilyFixed(tpms_family);
    // Disable contacts within drum components
    DEMSim.DisableContactBetweenFamilies(tpms_family, tpms_family);

////////////////////////////////////////////////////

    DEMSim.SetFamilyPrescribedPosition(recycled_family, "none", "Y+12", "none");
    DEMSim.SetFamilyPrescribedLinVel(recycled_family, "0", "none", "0");


    // DEMSim.SetFamilyPrescribedPosition(rec, "none", "(Y < -3) ? Y+9:Y", "none");
    // DEMSim.SetFamilyPrescribedLinVel(0, "(Y < -3) ? 0:vX", "none", "(Y < -3) ? 0:vZ");

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
    DEMSim.SetErrorOutVelocity(1e6);
    DEMSim.DisableJitifyClumpTemplates();

    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/DemoOutput_phx_tpms_periodic_10sec";
    create_directory(out_dir);

    float time_end = 10.0;
    unsigned int fps = 50;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    std::pair<float, float> plate_x_range = std::make_pair(-bxDim / 2., bxDim / 2.);
    std::pair<float, float> plate_z_range = std::make_pair(-bzDim / 2., bzDim / 2.);
    std::pair<float, float> plate_y_range = std::make_pair(-6., -3.);

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

            int num_changed = DEMSim.ChangeClumpFamily(recycled_family, plate_x_range, plate_y_range, plate_z_range);
            std::cout << "number of family changed: " << num_changed << std::endl;
            DEMSim.DoDynamicsThenSync(0.);



        }

        DEMSim.DoDynamics(step_size);

        if (curr_step % out_steps == 0){
            DEMSim.ChangeFamily(recycled_family, 0);
            DEMSim.DoDynamicsThenSync(0.);

        }


    }

    char cp_filename[200];
    sprintf(cp_filename, "%s/more_particles_settled_3e5.csv", out_dir.c_str());
    DEMSim.WriteClumpFile(std::string(cp_filename));

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    char cnt_filename[200];
    sprintf(cnt_filename, "%s/Contact_pairs_3e5.csv", out_dir.c_str());
    DEMSim.WriteContactFile(std::string(cnt_filename));


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count()) / time_end * 10.0 << " seconds (wall time) to finish 3 seconds' simulation"
              << std::endl;
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    std::cout << "DEMdemo_Centrifuge exiting..." << std::endl;
    return 0;
}
