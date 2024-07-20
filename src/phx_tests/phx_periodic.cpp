//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// 
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

// define a union of pin types, cylinder, half_teardop, full_teardrop
enum PinType {
    CYLINDER = 0,
    HALF_TEARDROP = 1,
    FULL_TEARDROP = 2
};

// Read pin positions from a file
// std::vector<float3> ReadPinPositions(std::string filename) {
//     std::vector<float3> pin_positions;

//     std::ifstream pin_pos_file(GetDEMEDataFile(filename));
//     std::string line;
//     std::getline(pin_pos_file, line); // skip the first line
//     while (std::getline(pin_pos_file, line)) {
//         std::istringstream ss(line);
//         std::string token;
//         std::vector<std::string> tokens;
//         while (std::getline(ss, token, ',')) {
//             tokens.push_back(token);
//         }
//         float3 pin_center = make_float3(std::stof(tokens[0]), std::stof(tokens[1]), std::stof(tokens[2]));
//         pin_positions.push_back(pin_center);

//     }
        
//     return pin_positions;
// }


int main(int argc, char* argv[]) {

    PinType pin_type = PinType::CYLINDER;
    std::string input_particle_positions;

    switch (pin_type){
        case PinType::CYLINDER:
            input_particle_positions = "../../input_particle_positions/cylinder_pins_particles_settled.csv";
            break;
        
        case PinType::HALF_TEARDROP:
            input_particle_positions = "../../input_particle_positions/teardrop_half_particles_settled.csv";
            break;

        case PinType::FULL_TEARDROP:
            input_particle_positions = "../../input_particle_positions/teardrop_full_particles_settled.csv";
            break;

        default:
            break;
    }

    float coeff_res = 0.7;  // from https://doi.org/10.1063/5.0149000

    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;


    int plate_family = 2;
    int sand_family = 0;
    int recylcled_family = 1;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL | OUTPUT_CONTENT::FAMILY);
    // DEMSim.SetVerbosity(STEP_METRIC);

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    DEMSim.SetNoForceRecord();

    auto mat_type_sand = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", coeff_res}, {"mu", 0.5}, {"Crr", 0.01}});
    auto mat_type_plate = DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.3}, {"CoR", coeff_res}, {"mu", 0.5}, {"Crr", 0.01}});


    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.}, 
                                      {-36.       , byDim},
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
    double sand_density = 3.6e3;
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

    ///////////////////////////////////////////////////////////////////////
    // Cylindrical pins //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    // read pin center from data file in data/sim_data/pin_pos.csv, skip the first header row, and assemble the value into a list of float3
    std::vector<float3> pin_centers = ReadPinPositions("sim_data/pin_pos.csv");
    std::cout << "number of pins: " << pin_centers.size() << std::endl;

    for (auto pin_center : pin_centers) {

        switch (pin_type){
            case PinType::CYLINDER:
{                float3 CylAxis = make_float3(0, 0, 1);
                float tube_radius = 0.2;
                auto pin = DEMSim.AddExternalObject();
                pin->AddCylinder(pin_center, CylAxis, tube_radius, mat_type_plate, 1); // outward normal

                break;
}                
            case PinType::HALF_TEARDROP:
{                auto pin = DEMSim.AddWavefrontMeshObject("../data/mesh/half_teardrop.obj", mat_type_plate);
                float4 rot = make_float4(1, 0, 0, 0);
                pin->Move(pin_center, rot);
                break;
}
            case PinType::FULL_TEARDROP:
{                auto pin = DEMSim.AddWavefrontMeshObject("../data/mesh/full_teardrop.obj", mat_type_plate);
                float4 rot = make_float4(1, 0, 0, 0);
                pin->Move(pin_center, rot);
                break;
}
            default:
                break;
        }
        std::cout << "added pin at " << pin_center.x << ", " << pin_center.y << ", " << pin_center.z << std::endl;

    }

    // if (pin_type == PinType::HALF_TEARDROP || pin_type == PinType::FULL_TEARDROP){
    //     DEMSim.SetFamilyFixed(10);
    // }

    // load the particles

    auto part1_clump_xyz = DEMSim.ReadClumpXyzFromCsv(input_particle_positions);
    auto part1_clump_quaternion = DEMSim.ReadClumpQuatFromCsv(input_particle_positions);

    std::vector<float3> in_xyz;
    std::vector<float4> in_quat;
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    unsigned int t_num = 0;
    unsigned int sphere_types = 4;
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
    base_batch.SetFamily(sand_family);
    DEMSim.AddClumps(base_batch);


    ///////////////// Load bottom plate particles /////////////////////////

    // Read plate particles from data file in data/sim_data/plate_pos.csv, skip the first header row, and assemble the value into a list of float3
    std::vector<float3> plate_particles;
    std::ifstream plate_pos_file(GetDEMEDataFile("clumps/bottom_plate_15mm_one_hole.csv"));
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
    std::cout << plate_particles.size() << " spheres make up the bottom plate" << std::endl;

    // Add drum
    auto plate = DEMSim.AddClumps(plate_template, make_float3(0));
    plate->SetFamilies(plate_family);

    DEMSim.SetFamilyFixed(plate_family);
    // Disable contacts within drum components
    DEMSim.DisableContactBetweenFamilies(plate_family, plate_family);

    // prescribe motion
    DEMSim.SetFamilyPrescribedPosition(recylcled_family, "none", "Y+26", "none");
    DEMSim.SetFamilyPrescribedLinVel(recylcled_family, "0", "none", "0");

////////////////////////////////////////////////////

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
    DEMSim.SetErrorOutVelocity(1e8);
    DEMSim.Initialize();

    path out_dir = current_path();

    out_dir += "/phx_cylinder_one_hole";

    create_directory(out_dir);

    float time_end = 10.;
    unsigned int fps = 100;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    unsigned int write_fps = 2000; // write 2000 frames per second 
    unsigned int write_out_steps = (unsigned int)(1.0 / (write_fps * step_size));
    unsigned int csv_frame = 0;

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();


    std::pair<float, float> plate_x_range = std::make_pair(-bxDim / 2., bxDim / 2.);
    std::pair<float, float> plate_z_range = std::make_pair(-bzDim / 2., bzDim / 2.);
    std::pair<float, float> plate_y_range = std::make_pair(-28., -26.);

    for (double t = 0; t < (double)time_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            DEMSim.ShowThreadCollaborationStats();
            max_v = max_v_finder->GetValue();
            std::cout << "Max velocity of any point in simulation is " << max_v << std::endl;

            //////////////////periodic boundary
            int num_changed = DEMSim.ChangeClumpFamily(recylcled_family, plate_x_range, plate_y_range, plate_z_range);
            std::cout << "number of family changed: " << num_changed << std::endl;
            DEMSim.DoDynamicsThenSync(0.);

            std::cout << "Frame: " << currframe << std::endl;
            std::cout << "Time: " << t << std::endl;

            char filename[200];
            sprintf(filename, "%s/DEM_frame_%04d.csv", out_dir.c_str(), csv_frame);
            DEMSim.WriteSphereFile(std::string(filename));
            std::cout << "write file: " << filename << std::endl;
            csv_frame++;
        }
        // only write the last second of data
        // if (curr_step % write_out_steps == 0 && t >= (double) time_end - 1 ){
        //     char filename[200];
        //     sprintf(filename, "%s/DEM_frame_%04d.csv", out_dir.c_str(), csv_frame);
        //     DEMSim.WriteSphereFile(std::string(filename));
        //     std::cout << "write file: " << filename << std::endl;
        //     csv_frame++;
        // }
        DEMSim.DoDynamics(step_size);
        DEMSim.ChangeFamily(recylcled_family, sand_family);
        DEMSim.DoDynamicsThenSync(0.);
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
