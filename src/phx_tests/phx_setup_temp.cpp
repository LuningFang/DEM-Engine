//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// TODO: at the end of the script, print out stats like bulk density, 
// Settling particles in the test rig, with prescribed pin temperature
// based on vertical position 
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

std::string force_model();

int main() {


    const float specific_heat = 9e6;

    // modify this later! start with some number first 
    double init_temp_cyl_celcius = 130.f;
    double init_temp_cyl_kelvin = init_temp_cyl_celcius + 273.15;

    double init_temp_sand_celcius = 40;
    double init_temp_sand_kelvin = init_temp_sand_celcius + 273.15;


    double bxDim = 5.0;
    double byDim = 48.0;
    double bzDim = 0.5;

    DEMSolver DEMSim;
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // Output family numbers (used to identify the centrifuging effect)
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL | OUTPUT_CONTENT::FAMILY | OUTPUT_CONTENT::GEO_WILDCARD);
    // DEMSim.SetVerbosity(STEP_METRIC);

    // If you don't need individual force information, then this option makes the solver run a bit faster.
    DEMSim.SetNoForceRecord();

    float fric_coef = 0.6;

    auto mat_type_sand = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", fric_coef}, {"Crr", 0.0}});
    auto mat_type_drum = DEMSim.LoadMaterial({{"E", 2e8}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", fric_coef}, {"Crr", 0.0}});


    DEMSim.InstructBoxDomainDimension({-bxDim / 2., bxDim / 2.}, 
                                      {-byDim / 2., byDim},
                                      {-bzDim / 2., bzDim/  2.});


    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_drum);


    // what does top open do? 
    // DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Since two types of materials have the same mu, this following call does not change the default mu for their
    // interaction, it's still 0.5.
    DEMSim.SetMaterialPropertyPair("mu", mat_type_sand, mat_type_drum, fric_coef);

    // Create some random clump templates for the filling materials
    // An array to store these generated clump templates
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;
    double sand_density = 3.6; // g/cm^3
    double scaling = 0.1;  // for testing, actual particle scale is 0.1

    double radius_array[3] = {0.212 * scaling, 0.2 * scaling, 0.178 * scaling};

    // Those following lines are needed. We must let the solver know that those var names are history variable etc.
    auto my_force_model = DEMSim.DefineContactForceModel(force_model());

    my_force_model->SetMustHaveMatProp({"E", "nu", "CoR", "mu", "Crr"});
    my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr"});
    my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});
    // wildcard values, temp for temperature, Q for heat flux
    my_force_model->SetPerGeometryWildcards({"Temp", "Q"});

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

    // read pin center from data file in data/sim_data/pin_pos.csv
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


    unsigned int num_clumps = input_xyz.size();
    std::cout << "number of particles: " << num_clumps << std::endl;

    // Casually select from generated clump types
    for (unsigned int i = 0; i < num_clumps; i++) {
        input_template_type.push_back(clump_types.at(i % clump_types.size()));
        // Every clump type that has a unique mass, gets a unique family number
        family_code.push_back((i % clump_types.size()) / 2);
    }
    auto particles = DEMSim.AddClumps(input_template_type, input_xyz);
    particles->SetFamilies(family_code);

    particles->AddGeometryWildcard("Q", std::vector<float>(num_clumps, 0));
    particles->AddGeometryWildcard("Temp", std::vector<float>(num_clumps, init_temp_sand_kelvin));



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
    DEMSim.SetFamilyExtraMargin(0, radius_array[0] * 0.2);

    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/phx_settling_temperature";
    create_directory(out_dir);

    float time_end = 3.0;
    unsigned int fps = 100;
    double frame_time = 1./double(fps);
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

            // now update temperature
            std::vector<float> Q_values = DEMSim.GetSphereWildcardValue(0, "Q", num_clumps);
            // get average value in std::vector

            std::vector<float> T_values = DEMSim.GetSphereWildcardValue(0, "Temp", num_clumps);
            // This is where I'm going to update T based on Q values
            for (int i = 0; i < num_clumps; i++) {
                T_values[i] += Q_values[i] * frame_time / (DEMSim.GetOwnerMass(i) * specific_heat);

            }
            DEMSim.SetSphereWildcardValue(0, "Temp", T_values);
            DEMSim.SetSphereWildcardValue(0, "Q", std::vector<float>(num_clumps, 0.f));

            // compute average temp
            double avg_temp = 0;
            for (int i = 0; i < num_clumps; i++) {
                avg_temp += T_values[i];
            }
            avg_temp /= num_clumps;
            std::cout << "average temperature: " << avg_temp << std::endl;
        }
        DEMSim.DoDynamics(step_size);
    }

    char cp_filename[200];
    sprintf(cp_filename, "%s/3_types_cylinder.csv", out_dir.c_str());
    DEMSim.WriteClumpFile(std::string(cp_filename));


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << (time_sec.count())  << " seconds (wall time) to finish " << time_end << " seconds simulation"
              << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

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
const float ks = 2e5;

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

    // this is where particle fluid particle heat exchange model should be added
    int curr_step = (int) (time / ts);
    if (curr_step % 2000 == 0) {
        double distances[11] = {-0.200, -0.160, -0.120, -0.080, -0.040, 0.000, 0.040, 0.080, 0.120, 0.16, 0.200};
        double volumes[11] = {0.150, 0.167, 0.182, 0.197, 0.217, 0.409, 0.077, 0.046, 0.027, 0.012, 0.000};

        double interval = 0.04;
        double start_distance = -0.2;

        // calculate index
        double ratio = overlapDepth / (ARadius + BRadius);
        int i = (int) ((ratio - start_distance) / interval);

        // Check if the distance is within bounds
        if (i > 0 && i < 11 - 1) {

            double t = (ratio - distances[i]) / (distances[i + 1] - distances[i]);
            double volume = volumes[i] * (1 - t) + volumes[i + 1] * t;
            // temperature
            double T_j = Temp_B[BGeo];
            double T_i = Temp_A[AGeo];
            double Q_ij = ks * volume *  (T_j - T_i);
            atomicAdd(Q_A + AGeo, Q_ij);
            atomicAdd(Q_B + BGeo, -Q_ij);
        }

    }

/////////////////////////////////////////////////
// Now we add an extra electrostatic force
////////////////////////////////////////////////

    )V0G0N";
    return model;
}