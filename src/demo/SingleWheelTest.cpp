//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// TODO: Add description
// =============================================================================

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

// #include <include/DEMdemo_WheelDPSimplified.h>

using namespace deme;

const double math_PI = 3.1415927;

// Load clump type, set clump material properties and load particle checkpoint file
void PrepareParticles(DEMSolver &DEMSim, std::shared_ptr<DEMMaterial> mat_type_wheel, std::string clump_file, double world_size_x, double world_size_y, double world_size_z);

/*
DEM Families:
0: Movable Terrain
1: Fixed Terrain
10: Wheel with angular velocity
20: Wheel with angular and forward velocity
*/

int main(int argc, char *argv[]) {
    // Process input data
    if (argc != 3)
    {
        std::cerr << "Usage: ./SingleWheelTest <slip> <batch_dir_name>" << std::endl;
        return 1;
    }
    // Import slip and batch directory from CLI arguments
    double slip = std::atof(argv[1]);
    std::string out_dir = std::string(argv[2]) + "_slip_" + std::to_string(slip);

    // Create output data structure:
    std::string rover_dir = out_dir + "/rover";
    std::string particles_dir = out_dir + "/particles";
    std::filesystem::create_directory(std::filesystem::path(out_dir));                 // create subdirectory with data for this specific run
    std::filesystem::create_directory(std::filesystem::path(rover_dir)); 
    std::filesystem::create_directory(std::filesystem::path(particles_dir));
    std::cout << "Output directory: " << out_dir << std::endl;

    // initialize parameters json file
    // temporary: hardcoded json format, will improve in future to use libraries
    // std::filesystem::path output_params_path = out_dir / "params.json";
    // std::ofstream output_params;
    // output_params.open(output_params_path);
    // output_params << "{\n\t\"slip\": " + std::to_string(slip) + "\n}";
    // output_params.close();

    // initialize output data CSV file
    std::string output_filename = out_dir + "/output.csv";
    std::ofstream output_datafile(output_filename);
    output_datafile << "t, f_x, f_y, f_z, d_c, v_x, v_max" << std::endl;

    // Define mesh and clump files
    // auto wheel_mesh_file = "/usr/local/share/chrono/data/robot/moonranger/obj/moonranger_wheel.obj";
    // auto clump_file = "/home/moonshot-chrono/sims/DEM-Engine/build_3/DemoOutput_GRCPrep_Part2/GRC_3e6.csv";
    // auto clump_file = "/home/moonshot-chrono/sims/DEM-Engine/build_3/DemoOutput_GRCPrep_Part2/GRC_3e5_Reduced_Footprint.csv";
    auto clump_file = "GRC_3e5.csv";
    auto wheel_mesh_file = GetDEMEDataFile("mesh/rover_wheels/moonranger_wheel.obj");

    // Simulation Settings
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetContactOutputContent({"OWNER", "FORCE", "POINT"});

    // Family 1 is for fixed ground which does not participate the force calc.
    DEMSim.SetFamilyFixed(1);
    DEMSim.DisableContactBetweenFamilies(1, 1);


    // If you don't need individual force information, then this option makes the solver run a bit faster.
    DEMSim.SetNoForceRecord();

    auto mat_type_wheel = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.00}});

    // `World'
    float G_mag = 9.81;
    float step_size = 1e-6;
    double world_size_x = 1.;
    double world_size_y = 0.3;
    // double world_size_y = 0.7;
    // double world_size_y = 1.;
    double world_size_z = 2.;

    // Define the wheel geometry
    float wheel_radius = 0.085; 
    float wheel_width = 0.06;
    // float wheel_weight = 100.;
    // float wheel_mass = wheel_weight / G_mag;
    float wheel_mass = 0.238; // kg
    float total_mass = 4.5; // kg
    float total_pressure = total_mass * G_mag; // N
    std::cout << "Total Pressure: " << total_pressure << "N" << std::endl;
    float added_pressure = (total_mass - wheel_mass) * G_mag; // N
    std::cout << "Added Pressure: " << added_pressure << "N" << std::endl;
    float wheel_IYY = wheel_mass * wheel_radius * wheel_radius / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_radius * wheel_radius + wheel_width * wheel_width);
    auto wheel =
        DEMSim.AddWavefrontMeshObject(wheel_mesh_file, mat_type_wheel);
    wheel->SetMass(total_mass);
    wheel->SetMOI(make_float3(wheel_IXX, wheel_IYY, wheel_IXX));
    // Give the wheel a family number so we can potentially add prescription
    wheel->SetFamily(10);
    // Track it
    auto wheel_tracker = DEMSim.Track(wheel);
    
    PrepareParticles(DEMSim, mat_type_wheel, clump_file, world_size_x, world_size_y, world_size_z);

    // Families' prescribed motions
    float w_r = 0.2;
    float v_ref = w_r * wheel_radius;

    double sim_end = 45.;
    // Note: this wheel is not `dictated' by our prescrption of motion because it can still fall onto the ground
    // (move freely linearly)
    DEMSim.SetFamilyPrescribedAngVel(10, "0", to_string_with_precision(w_r), "0", false);
    // An extra force (acceleration) is addedd to simulate the load that the wheel carries
    // DEMSim.AddFamilyPrescribedAcc(1, "none", "none", to_string_with_precision(-added_pressure / wheel_mass));
    // std::cout << "Prescribing Acceleration: " << -added_pressure / wheel_mass << "m/s^2" << std::endl;

    // `Real sim' family number
    DEMSim.SetFamilyPrescribedAngVel(20, "0", to_string_with_precision(w_r), "0", false);
    // Note: this wheel is not `dictated' by our prescrption of motion (hence the false argument), because it
    // can sink into the ground (move on Z dir); but its X and Y motions are explicitly controlled.
    // This one says when the experiment is going, the slip ratio is 0.5 (by our prescribing linear and angular vel)
    DEMSim.SetFamilyPrescribedLinVel(20, to_string_with_precision(v_ref * (1-slip)), "0", "none", false);
    // An extra force (acceleration) is addedd to simulate the load that the wheel carries
    // DEMSim.AddFamilyPrescribedAcc(2, "none", "none", to_string_with_precision(-added_pressure / wheel_mass));

    // Some inspectors
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
    auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    // Make ready for simulation
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -G_mag));
    // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
    // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
    // happen anyway and if it does, something already went wrong.
    DEMSim.SetMaxVelocity(50.);
    // Error out vel is used to force the simulation to abort when something goes wrong.
    DEMSim.SetErrorOutVelocity(60.);
    DEMSim.SetExpandSafetyMultiplier(1.1);
    //// TODO: Implement the better CDUpdateFreq adapt algorithm that overperforms the current one...
    DEMSim.SetCDUpdateFreq(30);
    // DEMSim.DisableAdaptiveUpdateFreq();
    DEMSim.Initialize();

    unsigned int fps = 10;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    unsigned int curr_step = 0;
    unsigned int currframe = 0;
    double frame_time = 1.0 / fps;
    unsigned int report_ps = 1000;
    unsigned int report_steps = (unsigned int)(1.0 / (report_ps * step_size));
    std::cout << "Output at " << fps << " FPS" << std::endl;

    // Put the wheel in place, then let the wheel sink in initially
    float max_z = max_z_finder->GetValue();
    wheel_tracker->SetPos(make_float3(-0.25, 0, max_z + 0.03 + wheel_radius));
    for (double t = 0; t < 0.5; t += frame_time) {
        char filename[200], meshname[200];
        std::cout << "Outputting frame: " << currframe << std::endl;
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", particles_dir.c_str(), currframe);
        sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", rover_dir.c_str(), currframe++);
        std::cout << "Writing sphere" << std::endl;
        DEMSim.WriteSphereFile(std::string(filename));
        std::cout << "Writing mesh" << std::endl;
        DEMSim.WriteMeshFile(std::string(meshname));
        std::cout << "Doing dynamics" << std::endl;
        DEMSim.DoDynamicsThenSync(frame_time);
    }

    // Switch wheel from free fall into DP test
    DEMSim.DoDynamicsThenSync(0);
    // You don't have to do this! I am just testing if sync-ing it twice breaks the system.
    DEMSim.DoDynamicsThenSync(0);
    DEMSim.ChangeFamily(10, 20);
    step_size *= 2.;
    DEMSim.UpdateStepSize(step_size);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // Active box domain
    float box_halfsize_x = wheel_radius * 1.2;
    float box_halfsize_y = wheel_width * 2.;

    for (double t = 0; t < sim_end; t += step_size, curr_step++) {

        if (curr_step % out_steps == 0) {
            // Update Active Box Domain
            // DEMSim.DoDynamicsThenSync(0.);
            // // Fix all particles
            // DEMSim.ChangeClumpFamily(1);
            // size_t num_changed = 0;
            // std::pair<float, float> Xrange =
            //     std::pair<float, float>(wheel_tracker->Pos().x - box_halfsize_x, wheel_tracker->Pos().x + box_halfsize_x);
            // std::pair<float, float> Yrange =
            //     std::pair<float, float>(wheel_tracker->Pos().y - box_halfsize_y, wheel_tracker->Pos().y + box_halfsize_y);
            // num_changed += DEMSim.ChangeClumpFamily(0, Xrange, Yrange);
            // std::cout << num_changed << " particles changed family number." << std::endl;

            // Write data
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << currframe << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", particles_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", rover_dir.c_str(), currframe++);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));
            DEMSim.ShowThreadCollaborationStats();
        }

        if (curr_step % report_steps == 0) {
            float3 forces = wheel_tracker->ContactAcc();
            forces *= total_mass;
            // output to terminal
            std::cout << "Time: " << t << std::endl;
            std::cout << "Force on wheel: " << forces.x << ", " << forces.y << ", " << forces.z << std::endl;
            std::cout << "Drawbar pull coeff: " << forces.x / total_pressure << std::endl;
            // std::cout << "Max system velocity: " << max_v_finder->GetValue() << std::endl;
            // output to file
            output_datafile << t << "," << forces.x << ", " << forces.y << ", " << forces.z << ", " << forces.x/total_pressure << ", " << wheel_tracker->Vel().z <<", " << max_v_finder->GetValue() << std::endl;
            output_datafile.flush();
        }

        DEMSim.DoDynamics(step_size);
    }

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    std::cout << "WheelDPSimpilified demo exiting..." << std::endl;
    return 0;
}


void PrepareParticles(DEMSolver &DEMSim, std::shared_ptr<DEMMaterial> mat_type_wheel, std::string clump_file, double world_size_x, double world_size_y, double world_size_z)
{
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.00}});
    // If you don't have this line, then mu between drum material and granular material will be the average of the
    // two.
    DEMSim.SetMaterialPropertyPair("mu", mat_type_wheel, mat_type_terrain, 0.8);

    DEMSim.InstructBoxDomainDimension(world_size_x, world_size_y, world_size_z);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);
    float bottom = -0.5;
    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
    // DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_wall);
    DEMSim.AddBCPlane(make_float3(0, world_size_y / 2, 0), make_float3(0, -1, 0), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0, -world_size_y / 2, 0), make_float3(0, 1, 0), mat_type_terrain);
    // X-dir bounding planes
    DEMSim.AddBCPlane(make_float3(-world_size_x / 2., 0, 0), make_float3(1, 0, 0), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(world_size_x / 2., 0, 0), make_float3(-1, 0, 0), mat_type_terrain);

    // Define the terrain particle templates
    // Calculate its mass and MOI
    float terrain_density = 2.6e3;
    float volume1 = 4.2520508;
    float mass1 = terrain_density * volume1;
    float3 MOI1 = make_float3(1.6850426, 1.6375114, 2.1187753) * terrain_density;
    float volume2 = 2.1670011;
    float mass2 = terrain_density * volume2;
    float3 MOI2 = make_float3(0.57402126, 0.60616378, 0.92890173) * terrain_density;
    // Scale the template we just created
    std::vector<double> scales = {0.0014, 0.00075833, 0.00044, 0.0003, 0.0002, 0.00018333, 0.00017};
    std::for_each(scales.begin(), scales.end(), [](double &r)
                  { r *= 10.; });
    // Then load it to system
    std::shared_ptr<DEMClumpTemplate> my_template2 =
        DEMSim.LoadClumpType(mass2, MOI2, GetDEMEDataFile("clumps/triangular_flat_6comp.csv"), mat_type_terrain);
    std::shared_ptr<DEMClumpTemplate> my_template1 =
        DEMSim.LoadClumpType(mass1, MOI1, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {my_template2,
                                                                                DEMSim.Duplicate(my_template2),
                                                                                my_template1,
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1)};
    // Now scale those templates
    for (int i = 0; i < scales.size(); i++)
    {
        std::shared_ptr<DEMClumpTemplate> &my_template = ground_particle_templates.at(i);
        // Note the mass and MOI are also scaled in the process, automatically. But if you are not happy with this, you
        // can always manually change mass and MOI afterwards.
        my_template->Scale(scales.at(i));
        // Give these templates names, 0000, 0001 etc.
        char t_name[20];
        sprintf(t_name, "%04d", i);
        my_template->AssignName(std::string(t_name));
    }

    // Now we load clump locations from a checkpointed file
    std::cout << "Making terrain..." << std::endl;
    std::unordered_map<std::string, std::vector<float3>> clump_xyz;
    std::unordered_map<std::string, std::vector<float4>> clump_quaternion;
    try {
        clump_xyz = DEMSim.ReadClumpXyzFromCsv(clump_file);
        clump_quaternion = DEMSim.ReadClumpQuatFromCsv(clump_file);
    } catch (...) {
        throw std::runtime_error("You will need to finish the GRCPrep demos first to obtain the checkpoint file GRC_3e6.csv, "
                                 "in order to run this demo. That is a 4m by 2m GRC-1 terrain patch that is around 15cm "
                                 "thick. If you don't have access to it, you can go to the forum "
                                 "(https://groups.google.com/g/projectchrono) to ask the authors for it.");
    }
    std::vector<float3> in_xyz;                              // particle position
    std::vector<float4> in_quat;                             // particle quaternion
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types; // particle clump type
    unsigned int t_num = 0;
    for (int i = 0; i < scales.size(); i++) {
        char t_name[20];
        sprintf(t_name, "%04d", t_num);

        auto this_type_xyz = clump_xyz[std::string(t_name)];
        auto this_type_quat = clump_quaternion[std::string(t_name)];

        size_t n_clump_this_type = this_type_xyz.size();
        std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type
                  << std::endl;
        // Prepare clump type identification vector for loading into the system (don't forget type 0 in
        // ground_particle_templates is the template for rover wheel)
        std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                 ground_particle_templates.at(t_num));

        // Add them to the big long vector
        in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
        in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
        in_types.insert(in_types.end(), this_type.begin(), this_type.end());
        std::cout << "Added clump type " << t_num << std::endl;
        // Our template names are 0000, 0001 etc.
        t_num++;
    }

    // Now, we don't need all particles loaded, remove ones are not flat
    std::vector<notStupidBool_t> elem_to_remove(in_xyz.size(), 0);

    for (size_t i = 0; i < in_xyz.size(); i++) {
        if (in_xyz.at(i).z > -0.44)
            elem_to_remove.at(i) = 1;
    }
    in_xyz.erase(std::remove_if(in_xyz.begin(), in_xyz.end(),
                                [&elem_to_remove, &in_xyz](const float3& i) {
                                    return elem_to_remove.at(&i - in_xyz.data());
                                }),
                    in_xyz.end());
    in_quat.erase(std::remove_if(in_quat.begin(), in_quat.end(),
                                    [&elem_to_remove, &in_quat](const float4& i) {
                                        return elem_to_remove.at(&i - in_quat.data());
                                    }),
                    in_quat.end());
    in_types.erase(std::remove_if(in_types.begin(), in_types.end(),
                                    [&elem_to_remove, &in_types](const auto& i) {
                                        return elem_to_remove.at(&i - in_types.data());
                                    }),
                    in_types.end());

    // Finally, load the info into this batch
    DEMClumpBatch base_batch(in_xyz.size());
    base_batch.SetTypes(in_types);
    base_batch.SetPos(in_xyz);
    base_batch.SetOriQ(in_quat);

    DEMSim.AddClumps(base_batch);
}