//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// A meshed ball hitting a granular bed under gravity. A collection of different
// ball densities and drop heights are tested against loosely packed material.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <random>

using namespace deme;
using namespace std::filesystem;

void print(int id, const std::vector<int>& container) {
    std::cout << id << ". ";
    for (const int x : container)
        std::cout << x << ' ';
    std::cout << '\n';
}

void runDEME(std::string dir_output, float friction, float massMultiplier);
std::string force_model();

int main() {


    std::string tag = "test";

    float conctact_friction = 0.1;  // takes the value
    float massMultiplier = 100.00;     // takes the value

    std::string out_dir = "/DemoOutput_Hookean_mu_0.1_real/" ;
    out_dir += "Test_" + tag + "/";

    std::cout << "Running case with friction: " << conctact_friction << ", and Mass multiplier: " << massMultiplier
              << std::endl;
    std::cout << "Dir out is " << out_dir << std::endl;

    runDEME(out_dir, conctact_friction, massMultiplier);

    return 0;
}

void runDEME(std::string dir_output, float frictionMaterial, float massMultiplier) {
    DEMSolver DEMSim;
    DEMSim.UseFrictionalHertzianModel();
    DEMSim.SetVerbosity("ERROR");
    DEMSim.SetOutputFormat("CSV");
    DEMSim.SetOutputContent({"ABSV"});
    DEMSim.SetMeshOutputFormat("VTK");
    DEMSim.SetContactOutputContent(DEME_POINT | OWNER | FORCE | CNT_WILDCARD);

    path out_dir = current_path();
    out_dir += dir_output;
    remove_all(out_dir);
    create_directories(out_dir);

    double kn_ratio = 4e6;
    float terrain_rad = 0.01;
    float gravityMagnitude = 10;
    float sphere_mass = 0.01;


    double kn = kn_ratio * sphere_mass * gravityMagnitude / terrain_rad;
    double kt = kn;

    auto mat_type_sphere = DEMSim.LoadMaterial({{"kn", kn}, {"kt", kt}, {"mu", frictionMaterial}, {"CoR", 0.01}});
    auto mat_type_box = DEMSim.LoadMaterial({{"kn", kn}, {"kt", kt}, {"mu", frictionMaterial}, {"CoR", 0.01}});
    auto my_force_model = DEMSim.DefineContactForceModel(force_model());

    my_force_model->SetMustHaveMatProp({"kn", "kt", "mu", "CoR"});
    my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});

    float DTc = PI * terrain_rad * sqrt(1e3 / (1e7 / (2 * (1 + 0.33)))) / (0.8766 + 0.163 * 0.33);

    float step_size = 1e-7;
    std::cout << "step size is : " << step_size << std::endl;
    double world_sizeX = 122.0 * terrain_rad;  // 122 works fine for frictionless
    double world_sizeZ = 27 * terrain_rad;

    DEMSim.InstructBoxDomainDimension({-world_sizeX / 2., world_sizeX / 2.}, {-5 * terrain_rad, 5 * terrain_rad},
                                      {-1 * world_sizeZ, 10 * terrain_rad});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_box);

    // creating the two clump templates we need, which are just spheres
    std::vector<std::shared_ptr<DEMClumpTemplate>> templates_terrain;

    templates_terrain.push_back(DEMSim.LoadSphereType(sphere_mass, terrain_rad, mat_type_sphere));
    templates_terrain.push_back(DEMSim.LoadSphereType(sphere_mass, terrain_rad, mat_type_sphere));

    // templates_terrain=DEMSim.LoadSphereType();
    unsigned int num_particle = 0;

    auto data_xyz = DEMSim.ReadClumpXyzFromCsv("../data/clumps/xyz.csv");
    std::vector<float3> input_xyz;

    std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;

    std::cout << data_xyz.size() << std::endl;
    for (unsigned int i = 0; i < data_xyz.size(); i++) {
        char t_name[20];
        sprintf(t_name, "%d", i);

        auto this_type_xyz = data_xyz[std::string(t_name)];
        input_xyz.insert(input_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());

        input_pile_template_type.push_back(templates_terrain[0]);
    }

    std::cout << input_xyz.size() << std::endl;
    std::cout << input_pile_template_type.size() << std::endl;

    auto allParticles = DEMSim.AddClumps(input_pile_template_type, input_xyz);
    allParticles->SetFamily(1);

    auto zeroParticle = DEMSim.AddClumps(templates_terrain[1], make_float3(0, 0, -terrain_rad));
    zeroParticle->SetFamily(3);
    auto driver = DEMSim.Track(zeroParticle);

    float Aext = -gravityMagnitude * (massMultiplier);
    float timeApplication = abs(Aext) > abs(gravityMagnitude) ? 2 * sqrt(terrain_rad * abs(Aext))
                                                              : 2 * sqrt(terrain_rad * abs(gravityMagnitude));
    std::string Aext_pattern =
        to_string_with_precision(Aext) + "*erf(t/sqrt(" + to_string_with_precision(timeApplication) + "))";
    std::cout << "applying this force law " << timeApplication << std::endl;
    DEMSim.AddFamilyPrescribedAcc(2, "none", "none", Aext_pattern);

    num_particle += input_xyz.size();

    std::cout << "Total num of particles: " << (int)DEMSim.GetNumClumps() << std::endl;

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -gravityMagnitude));

    DEMSim.Initialize();

    float sim_time = 5.0;
    float time_settling = 2.0;
    unsigned int fps = 10;
    float frame_time = 1.0 / fps;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    double terrain_max_z;

    bool status = true;
    bool changeMaterial = true;

    for (float t = 0; t < sim_time; t += frame_time) {
        std::cout << "Output file: " << currframe << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);

        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteContactFile(std::string(cnt_filename));
        currframe++;

        DEMSim.DoDynamicsThenSync(frame_time);

        if (t > time_settling / 2 && changeMaterial) {
            DEMSim.DoDynamicsThenSync(0);
            std::cout << "Including restitution coefficient" << std::endl;
            changeMaterial = false;
            DEMSim.SetFamilyClumpMaterial(1, mat_type_sphere);
        }

        if (t > time_settling && status) {
            DEMSim.DoDynamicsThenSync(0);
            DEMSim.ChangeFamily(3, 2);
            std::cout << "Extra mass applied" << std::endl;
            status = false;
            DEMSim.DoDynamicsThenSync(timeApplication);
        }
    }

    DEMSim.ShowTimingStats();

    std::cout << "==============================================================" << std::endl;

    std::cout << "DEMdemo_2DForce exiting..." << std::endl;
}

std::string force_model() {

    std::string model = R"V0G0N(

        if (overlapDepth > 0) {
            // Material properties
            float kn_cnt, kt_cnt, gn_cnt, gt_cnt;
            float mu_cnt, CoR_cnt;

            {
                // E and nu are associated with each material, so obtain them this way
                float kn_A = kn[bodyAMatType];
                float kt_A = kt[bodyAMatType];
                float mu_A = mu[bodyAMatType];
                float cor_A = CoR[bodyAMatType];

                // retrieve the material properties of body B
                float kn_B = kn[bodyBMatType];
                float kt_B = kt[bodyBMatType];
                float mu_B = mu[bodyBMatType];
                float cor_B = CoR[bodyBMatType];

                // pick the max 
                kn_cnt = (kn_A > kn_B) ? kn_A : kn_B;
                kt_cnt = (kt_A > kt_B) ? kt_A : kt_B;
                mu_cnt = (mu_A > mu_B) ? mu_A : mu_B;
                CoR_cnt = (cor_A > cor_B) ? cor_A : cor_B;


                // compute gn_cnt and gt_cnt based on cor, etc
                const float loge = (CoR_cnt < DEME_TINY_FLOAT) ? log(DEME_TINY_FLOAT) : log(CoR_cnt);

                float mass_eff = (AOwnerMass * BOwnerMass) / (AOwnerMass + BOwnerMass);

                float tmp_g = 1. + deme::PI_SQUARED / (loge * loge);
                gn_cnt = sqrt( 4 * mass_eff * kn_cnt / tmp_g);
                gt_cnt = gn_cnt;

            }

            float3 rotVelCPA, rotVelCPB;
            {
                // We also need the relative velocity between A and B in global frame to use in the damping terms
                // To get that, we need contact points' rotational velocity in GLOBAL frame
                // This is local rotational velocity (the portion of linear vel contributed by rotation)
                rotVelCPA = cross(ARotVel, locCPA);
                rotVelCPB = cross(BRotVel, locCPB);
                // This is mapping from local rotational velocity to global
                applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA.x, rotVelCPA.y, rotVelCPA.z, AOriQ.w, AOriQ.x, AOriQ.y, AOriQ.z);
                applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB.x, rotVelCPB.y, rotVelCPB.z, BOriQ.w, BOriQ.x, BOriQ.y, BOriQ.z);
            }

            // A few re-usables
            float3 vrel_tan;
            float3 delta_tan = make_float3(delta_tan_x, delta_tan_y, delta_tan_z);

            // Normal force part
            {
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

                // normal force
                force += (kn_cnt * overlapDepth - gn_cnt * projection) * B2A;
                // printf("normal force: %f, %f, %f\n", force.x, force.y, force.z);
            }

            // Tangential force part
            if (mu_cnt > 0.0) {
                float3 tangent_force = kt_cnt * delta_tan + gt_cnt * vrel_tan;
                const float ft = length(tangent_force);
                if (ft > DEME_TINY_FLOAT) {
                    // Reverse-engineer to get tangential displacement
                    const float ft_max = length(force) * mu_cnt;
                    if (ft > ft_max) {
                        tangent_force = (ft_max / ft) * tangent_force;
                        delta_tan = (tangent_force - gt_cnt * vrel_tan) / (kt_cnt);
                    }
                } else {
                    tangent_force = make_float3(0, 0, 0);
                }
                // Use force to collect tangent_force
                force -= tangent_force;
                // printf("tangent force: %f, %f, %f\n", tangent_force.x, tangent_force.y, tangent_force.z);
            }

            // Finally, make sure we update those wildcards (in this case, contact history)
            delta_tan_x = delta_tan.x;
            delta_tan_y = delta_tan.y;
            delta_tan_z = delta_tan.z;
        } else { 
            // This is to be more rigorous. If in fact no physical contact, then contact wildcards (such as contact history)
            // should be cleared (they will also be automatically cleared if the contact is no longer detected).
            delta_time = 0;
            delta_tan_x = 0;
            delta_tan_y = 0;
            delta_tan_z = 0;
        }
)V0G0N";

    return model;
}