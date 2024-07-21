//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Carbo hsp DEM properties
// from https://doi.org/10.1063/5.0149000
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

int main(int argc, char* argv[]) {

    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;

    if (argc != 2) {
        std::cout << "Usage: ./phx_periodic <TestID> <orifice size (mm)>" << std::endl;
        return 1;
    }

    int Test_ID = std::stoi(argv[1]);
    int plate_family = 2;
    int sand_family = 0;
    int recylcled_family = 1;
    int pin_family = 10;

    std::string out_dir = "July_Run_" + std::to_string(Test_ID) + "/";

    std::string input_particle_positions = out_dir + "settling/settled.csv";
    std::string pin_pos_file = "sim_data/pin_pos_run" + std::to_string(Test_ID) + ".csv";
 
    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL | OUTPUT_CONTENT::FAMILY);

    auto mat_type_carbo = DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.6}, {"Crr", 0.0}});
    auto mat_type_wall = DEMSim.LoadMaterial({{"E", 2e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.6}, {"Crr", 0.0}});


    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.}, 
                                      {-36.       , byDim},
                                      {-bzDim / 2., bzDim/  2.});


    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_wall);

    double carbo_density = 3.6;  // g/cm^3
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    std::vector<double> radius_array = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};
    LoadParticlesFromFile(DEMSim, input_particle_positions, mat_type_carbo, radius_array, carbo_density, sand_family);

    ///////////////// Add pins /////////////////////
    // mesh pins
    if (Test_ID != 0){
        std::string pin_mesh_filename = "mesh/pins/teardrop_run" + std::to_string(Test_ID) + ".obj";
        AddMeshPins(DEMSim, pin_pos_file, pin_mesh_filename, mat_type_wall);
    }
    else 
        AddCylindricalPins(DEMSim, pin_pos_file, 0.2, mat_type_wall);

    ///////////////// Add orifice /////////////////////
    std::string orifice_filename = "clumps/bottom_plate_8mm.csv";
    AddOrificeParticles(DEMSim, orifice_filename, mat_type_wall, plate_family);

    // set up periodic boundary, where particles discharged added to the top
    DEMSim.SetFamilyPrescribedPosition(recylcled_family, "none", "Y+26", "none");
    DEMSim.SetFamilyPrescribedLinVel(recylcled_family, "0", "none", "0");

    float step_size = 2e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, -981, 0));
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.SetExpandSafetyAdder(6.0);
    DEMSim.SetErrorOutVelocity(1e8);
    DEMSim.Initialize();

    out_dir += "/discharge/";

    create_directory(out_dir);

    float time_end = 3.;
    unsigned int fps = 100;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
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
            
            // change family of particles that are discharged
            int num_changed = DEMSim.ChangeClumpFamily(recylcled_family, plate_x_range, plate_y_range, plate_z_range);
            std::cout << "number of family changed: " << num_changed << std::endl;
            DEMSim.DoDynamicsThenSync(0.);

            std::cout << "Frame: " << currframe << std::endl;
            std::cout << "Time: " << t << std::endl;
            float recycled_family_mass = DEMSim.GetFamilyMass(recylcled_family);
            std::cout << "Recycled family mass: " << recycled_family_mass << " g, mass flow rate: " <<  recycled_family_mass / (1./fps) << "g / sec"   << std::endl;

            char filename[200];
            sprintf(filename, "%s/DEM_frame_%04d.csv", out_dir.c_str(), csv_frame);
            DEMSim.WriteSphereFile(std::string(filename));
            std::cout << "write file: " << filename << std::endl;
            csv_frame++;
        }

        DEMSim.DoDynamics(step_size);
        DEMSim.ChangeFamily(recylcled_family, sand_family);
        DEMSim.DoDynamicsThenSync(0.);
    }

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count()) << " seconds (wall time) to finish " << time_end << " seconds' simulation" << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();
    return 0;
}
