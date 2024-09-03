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


double init_temp_array[5] = {18.5, 19, 19.3, 21.3, 21.1};
double backplate_temp_intercept_array[5] = {127.18, 80.038, 81.843, 78.379, 78.413};
double backplate_temp_slope_array[5] = {-3.4086, -1.3823, -1.1256, -0.9396, -0.7764};
double specific_heat_array[5] = {7.917e6, 7.353e6, 7.329e6, 7.329e6, 7.32e6};
std::string orifice_filename_array[5] = {"175e-2mm", "4mm", "6mm", "8mm", "10mm"};
// Model that describes the temperature of the system
std::string force_model(double Q_fpf_ratio = 1.0, double backplate_temp_slope = 0.0, double backplate_temp_intercept = 80.0);

int main(int argc, char* argv[]) {

    if (argc != 3){
        std::cout << "Usage: ./phx_temp_validation <TestID 1-5> <Q_fpf_ratio>" << std::endl;
        return 1;
    }

    int TestID = std::stoi(argv[1]);
    double Q_fpf_ratio = std::stod(argv[2]);

    double specific_heat = specific_heat_array[TestID-1];
    double init_temp_sand = init_temp_array[TestID-1];
    double backplate_temp_slope = backplate_temp_slope_array[TestID-1];
    double backplate_temp_intercept = backplate_temp_intercept_array[TestID-1];
    std::string orifice_filename = "clumps/validation_bottom_plate_" + orifice_filename_array[TestID-1] + ".csv";

    bool initialize_particle_temp = false; // do i initialize particle temperature in the beginning?

    // const float specific_heat = 7.917e6;
    // double init_temp_sand = 18.5;
    // double backplate_temp_slope = -3.41;
    // double backplate_temp_intercept = 127.18;
    // std::string orifice_filename = "clumps/validation_bottom_plate_15e-1mm.csv";


    double init_temp_cyl = 134.7;
    std::string out_dir = "Aug_validation/";
    std::string input_particle_positions = out_dir + "settling/settled_parallel_plates.csv";

    // Append the formatted parameters to out_dir
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(1) << "no_pins_Test_" << TestID << "_Q_" << Q_fpf_ratio;
    out_dir += oss.str();

    // create directory
    std::filesystem::create_directories(out_dir);

    // first and foremost, wrtie the force_model string to a file named temperature_model.txt
    std::ofstream force_model_file(out_dir + "/temperature_model.txt");
    // write whatever in force_model() to the file
    force_model_file << force_model(Q_fpf_ratio, backplate_temp_slope, backplate_temp_intercept);
    force_model_file.close();

    // first row is time,mass_flow_rate,avg_outlet_temp
    std::ofstream info_file(out_dir + "/info.csv");
    info_file << "time,mass_flow_rate,avg_outlet_temp" << std::endl;

    float fric_coeff = 0.6;
    float coeff_res = 0.6;
    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;
    float3 CylAxis = make_float3(0, 0, 1);
    float tube_radius = 0.2;
    double carbo_density = 3.6;  // g/cm^3
    double scaling = 0.1;  // for testing, actual particle scale is 0.1
    float step_size = 2e-6;  // actual 2e-6


    int plate_family = 2;
    int sand_family = 0;
    int recylcled_family = 1;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL | OUTPUT_CONTENT::FAMILY | OUTPUT_CONTENT::GEO_WILDCARD);
    DEMSim.SetNoForceRecord();

    auto mat_type_carbo = DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.6}, {"Crr", 0.0}});
    auto mat_type_wall  = DEMSim.LoadMaterial({{"E", 2e7}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.6}, {"Crr", 0.0}});


    // DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.},
    //                                   {-36.       , byDim},
    //                                   {-bzDim / 2., bzDim/  2.});

    DEMSim.InstructBoxDomainDimension({-bxDim / 2. * 1.1, bxDim / 2. * 1.1},
                                      {-36. * 1.1       , byDim * 1.1},
                                      {-bzDim / 2. * 1.1, bzDim/  2. * 1.1});

    DEMSim.InstructBoxDomainBoundingBC("none", mat_type_wall);


    auto walls = DEMSim.AddExternalObject();
    // Add four walls
    walls->AddPlane(make_float3( bxDim/2, 0, 0), make_float3(-1, 0, 0), mat_type_wall);
    walls->AddPlane(make_float3(-bxDim/2, 0, 0), make_float3( 1, 0, 0), mat_type_wall);
    walls->AddPlane(make_float3(0, 0,  bzDim/2), make_float3(0, 0, -1), mat_type_wall);
    walls->AddPlane(make_float3(0, 0, -bzDim/2), make_float3(0, 0,  1), mat_type_wall);


    auto wall_tracker = DEMSim.Track(walls);

    auto my_force_model = DEMSim.DefineContactForceModel(force_model(Q_fpf_ratio, backplate_temp_slope, backplate_temp_intercept));

    // Those following lines are needed. We must let the solver know that those var names are history variable etc.
    my_force_model->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr"});
    my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr"});
    my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});
    // wildcard values, temp for temperature, Q for heat flux
    my_force_model->SetPerGeometryWildcards({"Temp", "Q"});

    // Add Particles
    std::vector<double> radius_array = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};
    auto particles = LoadParticlesFromFile(DEMSim, input_particle_positions, mat_type_carbo, radius_array, carbo_density, sand_family);

    // Finally, load them into the system
    int num_particles = particles->GetNumClumps();
    std::cout << "number of particles added: " << num_particles << std::endl;
    particles->AddGeometryWildcard("Q", std::vector<float>(num_particles, 0));
    particles->AddGeometryWildcard("Temp", std::vector<float>(num_particles, init_temp_sand));
    


    ///////////////// Add orifice /////////////////////
    auto orifice_particles = AddOrificeParticles(DEMSim, orifice_filename, mat_type_wall, plate_family);
    // int num_orifice_particles = 904; // for 1.5mm orifice 952 particles
    // for 2mm orifice, 904 particles
    // int num_orifice_particles = 760; // for 5mm orifice 760 particles

    int num_orifice_particles = GetNumParticlesInFile(orifice_filename);
    // int num_orifice_particles = 616; // for 8mm orifice, 616 particles

    std::cout << "number of orifice particles added: " << num_orifice_particles << std::endl;
    orifice_particles->AddGeometryWildcard("Q", std::vector<float>(num_orifice_particles, 0));
    orifice_particles->AddGeometryWildcard("Temp", std::vector<float>(num_orifice_particles, init_temp_sand));


    // prescribe motion
    DEMSim.SetFamilyPrescribedPosition(recylcled_family, "none", "Y+22", "none");
    DEMSim.SetFamilyPrescribedLinVel(recylcled_family, "0", "none", "0");

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, -981, 0));
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.SetExpandSafetyAdder(6.0);
    DEMSim.SetErrorOutVelocity(1e8);
    DEMSim.SetFamilyExtraMargin(0, radius_array[0] * 0.2);  // Necessary for pfp heat transfer model
    DEMSim.Initialize();

    wall_tracker->SetGeometryWildcardValues("Temp", std::vector<float>(4, init_temp_cyl));
    wall_tracker->SetGeometryWildcardValues("Q", std::vector<float>(4, 0));



    float time_end = 20.;
    unsigned int fps = 100;
    double frame_time = 1./double(fps);
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;
    unsigned int csv_frame = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    std::pair<float, float> plate_x_range = std::make_pair(-bxDim / 2., bxDim / 2.);
    std::pair<float, float> plate_z_range = std::make_pair(-bzDim / 2., bzDim / 2.);
    std::pair<float, float> plate_y_range = std::make_pair(-23., -22.);


    // now we intiialize sand particle temperature based on their positions 

    if (initialize_particle_temp){
        std::vector<float> T_values = DEMSim.GetSphereWildcardValue(0, "Temp", num_particles + num_orifice_particles);
        // This is where I'm going to update T based owner position
        for (int i = 0; i < num_particles; i++) {
            T_values[i] = init_temp_sand + DEMSim.GetOwnerPosition(i).y * backplate_temp_slope;
        }
        DEMSim.SetSphereWildcardValue(0, "Temp", T_values);
    }


    for (double t = 0; t < (double)time_end; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            // DEMSim.ShowThreadCollaborationStats();
            //////////////////periodic boundary
            int num_changed = DEMSim.ChangeClumpFamily(recylcled_family, plate_x_range, plate_y_range, plate_z_range);
            DEMSim.DoDynamicsThenSync(0.);

            std::cout << "Time: " << t << std::endl;
            float recycled_family_mass = DEMSim.GetFamilyMass(recylcled_family);
            std::cout << "mass_flow_rate, " <<  recycled_family_mass / (1./fps)  << std::endl;
	    
            // write to info.csv
            info_file << t << "," << recycled_family_mass / (1./fps) << ",";
        }
	    if (curr_step % (out_steps * 10) == 0) {
            	char filename[200];
            	sprintf(filename, "%s/DEM_frame_%05d.csv", out_dir.c_str(), csv_frame);
            	DEMSim.WriteSphereFile(std::string(filename));
            	csv_frame++;
	    }

        DEMSim.DoDynamics(step_size);
        if (curr_step % out_steps == 0){
            DEMSim.ChangeFamily(recylcled_family, sand_family);
            DEMSim.DoDynamicsThenSync(0.);


            // now let's update Q and temp
            std::vector<float> Q_values = DEMSim.GetSphereWildcardValue(0, "Q", num_particles + num_orifice_particles);
            // get average value in std::vector

            std::vector<float> T_values = DEMSim.GetSphereWildcardValue(0, "Temp", num_particles + num_orifice_particles);
            int counters = 0;
            double avg_temp_outlet = 0;
            // This is where I'm going to update T based on Q values
            for (int i = 0; i < num_particles; i++) {
                T_values[i] += Q_values[i] * frame_time / (DEMSim.GetOwnerMass(i) * specific_heat);

                if (std::abs (DEMSim.GetOwnerPosition(i).y + 20.5) < 0.1)  {
                    counters ++;
                    avg_temp_outlet += T_values[i];
                }
                else if (DEMSim.GetOwnerPosition(i).y > -3 ) {
                    T_values[i] = init_temp_sand;                
                }
            }

            std::cout << "avg_outlet_temp, " << avg_temp_outlet / counters << std::endl;
            info_file << avg_temp_outlet / counters << std::endl;
            DEMSim.SetSphereWildcardValue(0, "Temp", T_values);
            DEMSim.SetSphereWildcardValue(0, "Q", std::vector<float>(num_particles + num_orifice_particles, 0.f));

            DEMSim.SetAnalWildcardValue(0, "Temp", std::vector<float>(4, init_temp_cyl));
            DEMSim.SetAnalWildcardValue(0, "Q", std::vector<float>(4, 0.f));
        }
    }

    info_file.close();
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count()) << " seconds (wall time) to finish " << time_end << "  seconds' simulation" << std::endl;
    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

    return 0;
}

std::string force_model(double Q_fpf_ratio, double backplate_temp_slope, double backplate_temp_intercept) {
    std::string model = R"V0G0N(
 /////////////////////////////////////////////////////////////
// The first part is just the standard full Hertzian--Mindlin
/////////////////////////////////////////////////////////////

// No need to do any contact force calculation if no contact. And it can happen,
// since we added extra contact margin for adding electrostatic forces before
// physical contacts emerge.
float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;
const float ks = 2e5; // conductivity of carbo particles
double E_real_over_dem = 3.4e3;


double backplate_temp_slope = )V0G0N" + std::to_string(backplate_temp_slope) + R"V0G0N(;
double backplate_temp_intercept = )V0G0N" + std::to_string(backplate_temp_intercept) + R"V0G0N(;
double Q_fpf_ratio = )V0G0N" + std::to_string(Q_fpf_ratio) + R"V0G0N(;

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



    // force magnitude
    double force_mag =  k_n * overlapDepth + gamma_n * projection;

    if (force_mag < 0) {
        force_mag = 0;
    }

    int curr_step = (int)(time / ts);

    // bottom plate owner family is 2, Q = 0 
    // Geo of front, and side walls are 0, 1, 2, no Q there
    if (curr_step % 5000 == 0 && AOwnerFamily != 2 && BOwnerFamily != 2 && AGeo != 0 && BGeo != 0 && AGeo != 1 && BGeo != 1 && AGeo != 2 && BGeo != 2) {

        // radius contact
        double radius_eff = (ARadius * BRadius) / (ARadius + BRadius);
        double E_heat_transfer = E_cnt * E_real_over_dem;
        double radius_contact = powf( 2.f * radius_eff / E_heat_transfer * force_mag, 1.0/3.0);

        // Initialize temperatures
        double T_i = 0;
        double T_j = 0;
        double Q_ij = 0;


        if (AGeo == 3 && (int)myContactType == 11) {
            // wall is bodyA, look up temperature of the wall based on particle B position
            T_i = backplate_temp_slope * BOwnerPos.y + backplate_temp_intercept;
            T_j = Temp_B[BGeo];
            Q_ij = 4. * ks * radius_contact * (T_j - T_i);
        } else if (BGeo == 3 && (int)myContactType == 11) {
            // wall is bodyB, look up temperature of the wall based on particle A position
            T_j = backplate_temp_slope * AOwnerPos.y + backplate_temp_intercept;
            T_i = Temp_A[AGeo];
            Q_ij = 4. * ks * radius_contact * (T_j - T_i);
        } else {
            // Both AGeo and BGeo are particles
            T_i = Temp_A[AGeo];
            T_j = Temp_B[BGeo];
            Q_ij = 2. * ks * radius_contact * (T_j - T_i);
        }
        
        atomicAdd(Q_A + AGeo,  Q_ij);
        atomicAdd(Q_B + BGeo, -Q_ij);

        if (T_i < 0 || T_j < 0) {
            printf("NEGATIVE TEMP! T_j: %f, T_i: %f, Q_ij: %f, radius_contact: %f, force_mag: %f, AGeo: %d, BGeo: %d\n", T_j, T_i, Q_ij, radius_contact, force_mag, AGeo, BGeo);
        }

    }




} else {
    // If no physical contact, then contact wildcards (such as contact history)
    // should be cleared (they will also be automatically cleared if the contact is no longer detected).
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;
    
}

    // particle fluid particle heat transfer model 
    // Note that this needs to happen all particles within the neighborhood
    int curr_step = (int) (time / ts);
    if (curr_step % 5000 == 0 && AOwnerFamily != 2 && BOwnerFamily != 2 && AGeo != 0 && BGeo != 0 && AGeo != 1 && BGeo != 1 && AGeo != 2 && BGeo != 2) {

        double distances[11] = {-0.200, -0.160, -0.120, -0.080, -0.040, 0.000, 0.040, 0.080, 0.120, 0.16, 0.200};
        double volumes[11] = {0.150, 0.167, 0.182, 0.197, 0.217, 0.409, 0.077, 0.046, 0.027, 0.012, 0.000};
        int num_pts = 11;
        double interval = 0.04;
        double start_distance = -0.2;

        // calculate index
        double ratio = overlapDepth / (ARadius + BRadius);
        int i = (int) ((ratio - start_distance) / interval);

        // Check if the distance is within bounds
        if (i > 0 && i < num_pts - 1) {

            double t = (ratio - distances[i]) / (distances[i + 1] - distances[i]);
            double volume = volumes[i] * (1 - t) + volumes[i + 1] * t;

            // Initialize temperatures
            double T_i = 0;
            double T_j = 0;
            double Q_ij = 0;

        if (AGeo == 3 && (int)myContactType == 11) {
            // wall is bodyA, look up temperature of the wall based on particle B position
            T_i = backplate_temp_slope * BOwnerPos.y + backplate_temp_intercept;
            T_j = Temp_B[BGeo];
            Q_ij = 2. * ks * volume *  (T_j - T_i) * Q_fpf_ratio;

        } else if (BGeo == 3 && (int)myContactType == 11) {
            // wall is bodyB, look up temperature of the wall based on particle A position
            T_j = backplate_temp_slope * AOwnerPos.y + backplate_temp_intercept;
            T_i = Temp_A[AGeo];
            Q_ij = 2. * ks * volume *  (T_j - T_i) * Q_fpf_ratio;
        } else {
            // Both AGeo and BGeo are particles
            T_i = Temp_A[AGeo];
            T_j = Temp_B[BGeo];
            Q_ij = ks * volume *  (T_j - T_i) * Q_fpf_ratio;

        }


            atomicAdd(Q_A + AGeo,  Q_ij);
            atomicAdd(Q_B + BGeo, -Q_ij);

            // note, I can't have negative temp, print all info when T is negative
            // if (T_i < 0 || T_j < 0) {
            //     printf("T_j: %f, T_i: %f, Q_ij: %f, volume: %f, AGeo: %d, BGeo: %d, ratio = %f, i = %d\n", T_j, T_i, Q_ij, volume, AGeo, BGeo, Q_fpf_ratio, i);
            // }

        }

    }

    )V0G0N";
    return model;
}