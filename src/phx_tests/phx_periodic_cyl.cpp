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


// Model that describes the temperature of the system
std::string force_model();

int main(int argc, char* argv[]) {
    const float specific_heat = 172.f;
    double init_temp_cyl = 278.f;
    double init_temp_sand = 348.f;


//    PinType pin_type = PinType::FULL_TEARDROP;
    std::string input_particle_positions = "../../input_particle_positions/cylinder_pins_particles_settled.csv";

    float fric_coeff = 0.6;
    float coeff_res = 0.6;
    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;


    int plate_family = 2;
    int sand_family = 0;
    int recylcled_family = 1;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL | OUTPUT_CONTENT::FAMILY | OUTPUT_CONTENT::GEO_WILDCARD);
    // DEMSim.SetVerbosity(STEP_METRIC);

    DEMSim.SetNoForceRecord();

    auto mat_type_sand = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", coeff_res}, {"mu", fric_coeff}, {"Crr", 0.01}});
    auto mat_type_plate = DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.3}, {"CoR", coeff_res}, {"mu", fric_coeff}, {"Crr", 0.01}});


    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.},
                                      {-36.       , byDim},
                                      {-bzDim / 2., bzDim/  2.});


    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_plate);

    DEMSim.SetMaterialPropertyPair("mu", mat_type_sand, mat_type_plate, fric_coeff);

    // Create some random clump templates for the filling materials
    // An array to store these generated clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    double sand_density = 3.6;
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    double radius_array[4] = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling, 0.125 * scaling};

    auto my_force_model = DEMSim.DefineContactForceModel(force_model());

    // Those following lines are needed. We must let the solver know that those var names are history variable etc.
    my_force_model->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr"});
    my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr"});
    my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});
    // wildcard values, temp for temperature, Q for heat flux
    my_force_model->SetPerGeometryWildcards({"Temp", "Q"});

    // ratio, 20%, 50%, 30%
    double mass;

    // Then randomly create some clump templates for filling the drum
    for (int i = 0; i < 4; i++) {
        double volume = 4./3. * PI * std::pow(radius_array[i],3);
        mass = volume * sand_density;
        // Load a sphere type
        clump_types.push_back(DEMSim.LoadSphereType(mass, radius_array[i], mat_type_sand));
    }

    // read pin center from data file in data/sim_data/pin_pos.csv, skip the first header row, and assemble the value into a list of float3
    std::vector<float3> pin_centers = ReadPinPositions("sim_data/pin_pos.csv");
    std::cout << "number of pins: " << pin_centers.size() << std::endl;


    auto pin = DEMSim.AddExternalObject();
    for (auto pin_center : pin_centers) {

        float3 CylAxis = make_float3(0, 0, 1);
        float tube_radius = 0.2;
        pin->AddCylinder(pin_center, CylAxis, tube_radius, mat_type_plate, 1); // outward normal
        std::cout << "added pin at " << pin_center.x << ", " << pin_center.y << ", " << pin_center.z << std::endl;
    }

    auto pin_tracker = DEMSim.Track(pin);

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
    auto particles = DEMSim.AddClumps(base_batch);
    int num_particles = in_xyz.size();
    particles->AddGeometryWildcard("Q", std::vector<float>(num_particles, 0));
    particles->AddGeometryWildcard("Temp", std::vector<float>(num_particles, init_temp_sand));




    ///////////////// Load bottom plate particles /////////////////////////

    // Read plate particles from data file in data/sim_data/plate_pos.csv, skip the first header row, and assemble the value into a list of float3
    std::vector<float3> plate_particles;
    std::ifstream plate_pos_file(GetDEMEDataFile("clumps/bottom_plate_8mm.csv"));
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

    int num_pins = pin_centers.size();
    pin_tracker->SetGeometryWildcardValues("Temp", std::vector<float>(num_pins, init_temp_cyl));
    pin_tracker->SetGeometryWildcardValues("Q", std::vector<float>(num_pins, 0));

    path out_dir = current_path();

    out_dir += "/phx_periodic_temp_cylinder";

    // add the coefficient of restitution to the output directory, use argv[1] as the folder name
    create_directory(out_dir);

    float time_end = 10.;
    unsigned int fps = 100;
    double frame_time = 1./double(fps);
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
        if (curr_step % out_steps == 0){
            DEMSim.ChangeFamily(recylcled_family, sand_family);
            DEMSim.DoDynamicsThenSync(0.);


            // now let's update Q and temp
            std::vector<float> Q_values = DEMSim.GetSphereWildcardValue(0, "Q", num_particles);
            // get average value in std::vector

            std::vector<float> T_values = DEMSim.GetSphereWildcardValue(0, "Temp", num_particles);
            // This is where I'm going to update T based on Q values
            for (int i = 0; i < num_particles; i++) {
                T_values[i] += Q_values[i] * frame_time / (DEMSim.GetOwnerMass(i) * specific_heat);
            }
            DEMSim.SetSphereWildcardValue(0, "Temp", T_values);
            DEMSim.SetSphereWildcardValue(0, "Q", std::vector<float>(num_particles, 0.f));

            // compute average temp
            double avg_temp = 0;
            for (int i = 0; i < num_particles; i++) {
                avg_temp += T_values[i];
            }
            avg_temp /= num_particles;
            std::cout << "average temperature: " << avg_temp << std::endl;

        }
    }

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();


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

std::string force_model() {
    std::string model = R"V0G0N(
 /////////////////////////////////////////////////////////////
// The first part is just the standard full Hertzian--Mindlin
/////////////////////////////////////////////////////////////

// No need to do any contact force calculation if no contact. And it can happen,
// since we added extra contact margin for adding electrostatic forces before
// physical contacts emerge.
float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;

if (overlapDepth > 0) {
    // Material properties
    // float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;
    {
        // E and nu are associated with each material, so obtain them this way
        float E_A = E[bodyAMatType];
        float nu_A = nu[bodyAMatType];
        float E_B = E[bodyBMatType];
        float nu_B = nu[bodyBMatType];
        matProxy2ContactParam<float>(E_cnt, G_cnt, E_A, nu_A, E_B, nu_B);
        // CoR, mu and Crr are pair-wise, so obtain them this way
        CoR_cnt = CoR[bodyAMatType][bodyBMatType];
        mu_cnt = mu[bodyAMatType][bodyBMatType];
        Crr_cnt = Crr[bodyAMatType][bodyBMatType];
    }

    float3 rotVelCPA, rotVelCPB;
        // We also need the relative velocity between A and B in global frame to use in the damping terms
        // To get that, we need contact points' rotational velocity in GLOBAL frame
        // This is local rotational velocity (the portion of linear vel contributed by rotation)
        rotVelCPA = cross(ARotVel, locCPA);
        rotVelCPB = cross(BRotVel, locCPB);
        // This is mapping from local rotational velocity to global
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA.x, rotVelCPA.y, rotVelCPA.z, AOriQ.w, AOriQ.x, AOriQ.y,
                                                AOriQ.z);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB.x, rotVelCPB.y, rotVelCPB.z, BOriQ.w, BOriQ.x, BOriQ.y,
                                                BOriQ.z);

    // A few re-usables
    float mass_eff, sqrt_Rd, beta;
    float3 vrel_tan;
    float3 delta_tan = make_float3(delta_tan_x, delta_tan_y, delta_tan_z);

    // Normal force part
        // The (total) relative linear velocity of A relative to B
        const float3 velB2A = (ALinVel + rotVelCPA) - (BLinVel + rotVelCPB);
        const float projection = dot(velB2A, B2A);
        vrel_tan = velB2A - projection * B2A;

        // Now we already have sufficient info to update contact history
        {
            delta_tan += ts * vrel_tan;
            const float disp_proj = dot(delta_tan, B2A);
            delta_tan -= disp_proj * B2A;
            delta_time += ts;
        }

        mass_eff = (AOwnerMass * BOwnerMass) / (AOwnerMass + BOwnerMass);
        sqrt_Rd = sqrt(overlapDepth * (ARadius * BRadius) / (ARadius + BRadius));
        const float Sn = 2. * E_cnt * sqrt_Rd;

        const float loge = (CoR_cnt < DEME_TINY_FLOAT) ? log(DEME_TINY_FLOAT) : log(CoR_cnt);
        beta = loge / sqrt(loge * loge + deme::PI_SQUARED);

        const float k_n = deme::TWO_OVER_THREE * Sn;
        const float gamma_n = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta * sqrt(Sn * mass_eff);

        force += (k_n * overlapDepth + gamma_n * projection) * B2A;

    // Rolling resistance part

    if (mu_cnt > 0.0) {
        const float kt = 8. * G_cnt * sqrt_Rd;
        const float gt = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta * sqrt(mass_eff * kt);
        float3 tangent_force = -kt * delta_tan - gt * vrel_tan;
        const float ft = length(tangent_force);
        if (ft > DEME_TINY_FLOAT) {
            // Reverse-engineer to get tangential displacement
            const float ft_max = length(force) * mu_cnt;
            if (ft > ft_max) {
                tangent_force = (ft_max / ft) * tangent_force;
                delta_tan = (tangent_force + gt * vrel_tan) / (-kt);
            }
        } else {
            tangent_force = make_float3(0, 0, 0);
        }
        // Use force to collect tangent_force
        force += tangent_force;
    }

    // Finally, make sure we update those wildcards (in this case, contact history)
    delta_tan_x = delta_tan.x;
    delta_tan_y = delta_tan.y;
    delta_tan_z = delta_tan.z;


    const float ks = 385.f;

    // force magnitude
    double force_mag =  k_n * overlapDepth + gamma_n * projection;

    if (force_mag < 0) {
        force_mag = 0;
    }

    int curr_step = (int)(time / ts);

    if (curr_step % 2000 == 0) {
        // radius contact
        double radius_eff = (ARadius * BRadius) / (ARadius + BRadius);
        double radius_contact = powf( 2.f * radius_eff / E_cnt * force_mag, 1.0/3.0);

        // temperature
        double T_j = Temp_B[BGeo];
        double T_i = Temp_A[AGeo];
        double Q_ij = 2. * ks * radius_contact *  (T_j - T_i);
        atomicAdd(Q_A + AGeo, Q_ij);
        atomicAdd(Q_B + BGeo, -Q_ij);

    }




} else {
    // This is to be more rigorous. If in fact no physical contact, then contact wildcards (such as contact history)
    // should be cleared (they will also be automatically cleared if the contact is no longer detected).
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;





}

/////////////////////////////////////////////////
// Now we add an extra electrostatic force
////////////////////////////////////////////////

    )V0G0N";
    return model;
}
