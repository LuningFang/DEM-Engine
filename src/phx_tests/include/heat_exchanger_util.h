// DEM_common_utils.h

#pragma once

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

using namespace deme;
using namespace std::filesystem;

// Read pin positions from a file
inline std::vector<float3> ReadPinPositions(const std::string& filename) {
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

// Read TPMS particle positions from a CSV file
std::vector<float3> ReadTPMSPositions(const std::string& filename) {
    std::vector<float3> positions;
    std::ifstream file(GetDEMEDataFile(filename));
    std::string line;
    // Skip the header
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::vector<std::string> result;
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            result.push_back(cell);
        }
        float3 pos = make_float3(std::stof(result[0]), std::stof(result[1]), std::stof(result[2]));
        positions.push_back(pos);
    }
    return positions;
}


// Populate particle positions given the domain size and the bounding box of the pin size 
inline std::vector<float3> PopulateParticlePositions(float spacing, float3 dim, const std::vector<float3>& pin_pos, double pin_radius) {
    std::vector<float3> material_points;
    PDSampler sampler(spacing);

    float box_X = dim.x;
    float box_Y = dim.y;
    float box_Z = dim.z;

    float fill_y_bottom = -box_Y / 2.0f + 2 * spacing;
    float fill_y_top = box_Y - spacing;
    float3 sampler_center = make_float3(0.0f, fill_y_bottom, 0.0f);
    float3 sampler_hdim = make_float3(box_X / 2.0f - spacing, 1e-6, box_Z / 2.0f - spacing);

    // find the location where the pin array starts and end
    std::vector<float> unique_y_pos;
    for (const auto& pin : pin_pos) {
        double y = pin.y;
        if (std::find(unique_y_pos.begin(), unique_y_pos.end(), y) == unique_y_pos.end()) {
            unique_y_pos.push_back(y);
        }
    }
    std::sort(unique_y_pos.begin(), unique_y_pos.end());
    int tube_pos_itr = 0;

    while (sampler_center.y <= fill_y_top) {
        if (tube_pos_itr >= unique_y_pos.size()) {
            break;
        }
        // get the upper and lower bound of the pins in each layer
        double curr_tube_lower_bound = unique_y_pos[tube_pos_itr] - pin_radius;
        double curr_tube_upper_bound = unique_y_pos[tube_pos_itr] + pin_radius;

        // now we add particles in between the pins
        if (sampler_center.y < curr_tube_lower_bound - 0.5 * spacing) {
            auto points = sampler.SampleBox(sampler_center, sampler_hdim);
            material_points.insert(material_points.end(), points.begin(), points.end());
            sampler_center.y += spacing;
        } else {
            sampler_center.y = curr_tube_upper_bound + 0.5 * spacing;
            tube_pos_itr++;
        }
    }

    while (sampler_center.y <= fill_y_top) {
        auto points = sampler.SampleBox(sampler_center, sampler_hdim);
        material_points.insert(material_points.end(), points.begin(), points.end());
        sampler_center.y += spacing;
    }

    return material_points;
}


// Add pins to the simulation (analytical circular pins)
std::vector<float3> AddCylindricalPins(DEMSolver& DEMSim, std::string pin_pos_filename, float pin_radius, std::shared_ptr<DEMMaterial> mat_type_wall) {
    std::vector<float3> pin_centers = ReadPinPositions(pin_pos_filename);
    std::cout << "number of pins: " << pin_centers.size() << std::endl;

    for (auto pin_center : pin_centers) {
        auto pin = DEMSim.AddExternalObject();
        pin->AddCylinder(pin_center, make_float3(0,0,1), pin_radius, mat_type_wall, 1); // outward normal
        std::cout << "added cylinder pin at " << pin_center.x << ", " << pin_center.y << ", " << pin_center.z << std::endl;
    }

    return pin_centers;
}

// Add pins to the simulation (mesh pins)
std::vector<float3> AddMeshPins(DEMSolver& DEMSim, std::string pin_pos_filename, std::string pin_obj_name, std::shared_ptr<DEMMaterial> mat_type_wall, int family_ID = 10) {
    std::vector<float3> pin_centers = ReadPinPositions(pin_pos_filename);
    std::cout << "number of pins: " << pin_centers.size() << std::endl;

    for (auto pin_center : pin_centers) {
        auto pin = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile(pin_obj_name), mat_type_wall);
        float4 rot = make_float4(0, 0, 0, 1);
        pin->Move(pin_center, rot);
        pin->SetFamily(family_ID);
        std::cout << "added meshed pin at " << pin_center.x << ", " << pin_center.y << ", " << pin_center.z << std::endl;
    }
    DEMSim.SetFamilyFixed(family_ID);
    return pin_centers;
}

std::shared_ptr<DEMClumpBatch> LoadParticlesFromFile(DEMSolver& DEMSim, 
                           const std::string& filename, std::shared_ptr<DEMMaterial> mat_type,
                           std::vector<double> radius_array,
                           float density,
                           int particle_family
                           ) {
    // An array to store clump templates based on radius and density
    std::vector<std::shared_ptr<DEMClumpTemplate>> clump_types;

    double mass;
    // Then randomly create some clump templates for filling the drum    
    for (int i = 0; i < 3; i++) {
        double volume = 4./3. * PI * std::pow(radius_array[i],3);
        mass = volume * density;
        // Load a sphere type
        clump_types.push_back(DEMSim.LoadSphereType(mass, radius_array[i], mat_type));
    }

    // load the particles
    auto clump_xyz = DEMSim.ReadClumpXyzFromCsv(filename);
    auto clump_quaternion = DEMSim.ReadClumpQuatFromCsv(filename);

    std::vector<float3> in_xyz;
    std::vector<float4> in_quat;
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    unsigned int t_num = 0;
    unsigned int sphere_types = radius_array.size();
    for (int i = 0; i < sphere_types; i++) {
        char t_name[20];
        sprintf(t_name, "%04d", t_num);

        auto this_type_xyz = clump_xyz[std::string(t_name)];
        auto this_type_quat = clump_quaternion[std::string(t_name)];

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
    base_batch.SetFamily(particle_family);
    auto particles = DEMSim.AddClumps(base_batch);
    return particles;
}

// Add bottom plate particles as a clump
std::shared_ptr<DEMClumpBatch> AddOrificeParticles(DEMSolver& DEMSim, 
                        const std::string& filename, std::shared_ptr<DEMMaterial> mat_type,
                        int plate_family = 20) {

    std::vector<float3> plate_particles;
    std::ifstream plate_pos_file(GetDEMEDataFile(filename));

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
                             std::vector<float>(plate_particles.size(), 0.025), plate_particles, mat_type);
    std::cout << plate_particles.size() << " spheres make up the bottom plate" << std::endl;

    // Add drum
    auto plate = DEMSim.AddClumps(plate_template, make_float3(0));
    plate->SetFamilies(plate_family);
    // Fix orifice particles
    DEMSim.SetFamilyFixed(plate_family);
    // Disable contacts within drum components
    DEMSim.DisableContactBetweenFamilies(plate_family, plate_family);

    return plate;

}