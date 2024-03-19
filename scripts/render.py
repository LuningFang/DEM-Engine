import bpy
import math
import random
import mathutils
import csv
import sys
import os.path
time = []


# Resolution of the output image
# res_x = 1920
res_x = 802
res_y = 1282

# look up table for clump type and corresponding radius
clump_radius = {'0000': 0.0212, '0001':0.02, '0002':0.0178, '0003': 0.025 }

print(sys.argv)
if len(sys.argv) != 5:
    print("usage ./blender --python --background animation_script_batch start_frame")
    sys.exit(1)

#============ specify the directory of csv, obj, image, and such
data_sim = "/srv/home/fang/phX/build_DEM/bin/DemoOutput_phx_periodic_8mm_orifice_2000fps/"

#data_sim = "C:/Users/fang/Documents/phx/dem_results/8mm_2000fps/"
image_dir = data_sim + "image/"

# if image_dir does not exist, create it
if not os.path.exists(image_dir):
    os.makedirs(image_dir)

#============ find the position that the camera points at
def point_at(obj, target, roll=0):
    """
    Rotate obj to look at target
    :arg obj: the object to be rotated. Usually the camera
    :arg target: the location (3-tuple or Vector) to be looked at
    :arg roll: The angle of rotation about the axis from obj to target in radians.
    """
    if not isinstance(target, mathutils.Vector):
        target = mathutils.Vector(target)
    loc = obj.location
    # direction points from the object to the target
    direction = target - loc

    quat = direction.to_track_quat('-Z', 'Y')

    # /usr/share/blender/scripts/addons/add_advanced_objects_menu/arrange_on_curve.py
    quat = quat.to_matrix().to_4x4()
    rollMatrix = mathutils.Matrix.Rotation(roll, 4, 'Z')

    # remember the current location, since assigning to obj.matrix_world changes it
    loc = loc.to_tuple()
    #obj.matrix_world = quat * rollMatrix
    # in blender 2.8 and above @ is used to multiply matrices
    # using * still works but results in unexpected behaviour!
    obj.matrix_world = quat @ rollMatrix
    obj.location = loc

#============ load position of each cylinder
dis = []

# read cylinder position from a file
cylinder_position_file = "cylinder_pos.csv"
# get x, y and z position, skip the first line, use csv to read the file
with open(cylinder_position_file, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    # skip the first row
    next(spamreader)
    for row in spamreader:
        dis.append((float(row[0]), float(row[1]), float(row[2])))

# total number of cylinders
num_cylinders = len(dis[0])
print(dis)


#===========================================
#============================ Start the loop
#===========================================
jobid = int(sys.argv[4])
start_frame = jobid * 100
end_frame = start_frame + 100
for i in range(start_frame, end_frame, 1):
    #===========================================
    #======== check if the png file exits or not
    #===========================================
    # image_path = image_dir + str(i) + ".png"
    # file_exists = os.path.exists(image_path)
    # if file_exists:
    #     sys.exit()

    #===========================================
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.objects.keys()

    file_loc = "cylinder_bld.obj"
    obj_name = "cylinder"
    obj_name_spe = "cylinder"

    for cylinder_pos in dis:

        bpy.ops.import_scene.obj(filepath=file_loc)
        # Get the newly imported object
        obj_object = bpy.context.selected_objects[0]

        # Check if an object was successfully imported
        if obj_object:
            print("cylinder pos: ", cylinder_pos)
            obj_object.name = obj_name  # Assign a unique name to each cylinder instance
            obj_object.location.x = cylinder_pos[0]
            obj_object.location.y = cylinder_pos[1]
            obj_object.location.z = cylinder_pos[2]

    # bpy.context.view_layer.update()

    bpy.context.view_layer.update()

    #===========================================
    #==================== Load SPH particle file
    #===========================================
    # load file name, padded with zeros
    particle_file_name = data_sim + "DEM_frame_{:06d}.csv".format(i)

    # read particle_file_name using pandas and skip the header
    # use csv to read the file, and append to positions array
    # make sure to skip the first line, which is the header
    count = 0
    count_blue = 0
    count_gray = 0
    positions_blue = []
    positions_gray = []
    radius_array = []
    for line in open(particle_file_name):
        if count == 0:
            count = count + 1
            continue
        else:
            # you have to parse "x", "y", "z" and "r" from the variable "line"
            line_seg = line.split(",")
            x, y, z = line_seg[0], line_seg[1], line_seg[2]

            clump_type = line_seg[7]

                    # if float(y) > 1.0 or float(y) < -1.0 or float(z) > 0.9:
            #     continue
            position_buff = (float(x), float(y), float(z))
            # color = line_seg[8]
            positions_gray.append(position_buff)
            radius_array.append(float(clump_radius[clump_type]))

            count_gray = count_gray + 1
            count = count + 1

    # print position of the last particle in blue
    print("count_gray: ", count_gray)
    print("particle position:")
    print(positions_gray[count_gray - 1])

    print("radius array size: " + str(len(radius_array)))
    """ -------------- PARTICLE SYSTEM START-------------- """
# Your positions_blue, positions_gray, count_blue, and radius_particle variables
# should be defined here.


# Your positions_blue, positions_gray, count_blue, and radius_particle variables
# should be defined here.

    context = bpy.context

    # Create the gray material
    material_gray = bpy.data.materials.new(name="GrayMaterial")
    material_gray.diffuse_color = (0.4, 0.4, 0.4, 1)  # Gray color

    # Create a separate icosphere for the gray particles and set the gray material
    bpy.ops.mesh.primitive_ico_sphere_add(radius=1, location=(50, 50, 50))
    ico_gray = context.object
    ico_gray.data.materials.append(material_gray)

    # Create a new cube for gray particles
    bpy.ops.mesh.primitive_cube_add(size=0.0001)
    gray_cube = context.object

    # Create the particle system for gray particles
    gray_ps = gray_cube.modifiers.new("ParticlesGray", 'PARTICLE_SYSTEM').particle_system
    gray_psname = gray_ps.name

    # Configure particle system settings for gray particles
    gray_ps.settings.count = len(positions_gray)
    gray_ps.settings.lifetime = 1000
    gray_ps.settings.frame_start = gray_ps.settings.frame_end = 1
    gray_ps.settings.render_type = "OBJECT"
    gray_ps.settings.instance_object = ico_gray  # Use the gray icosphere as the instance object

    # Clear the post frame handler
    bpy.app.handlers.frame_change_post.clear()

    # Register the particle_handler_gray with the post frame handler
    def particle_handler_gray(scene, depsgraph):
        ob = depsgraph.objects.get(gray_cube.name)
        if ob:
            gray_ps = ob.particle_systems[gray_psname]
            f = scene.frame_current
            for m, particle in enumerate(gray_ps.particles):
                setattr(particle, "location", positions_gray[m])
                setattr(particle, "size", radius_array[m])

    # Register both particle handlers
    bpy.app.handlers.frame_change_post.append(particle_handler_gray)

    # Trigger frame update
    bpy.context.scene.frame_current = 2
    bpy.context.view_layer.update()

    #===========================================
    #========================= Rendering setting
    #===========================================
    # bpy.ops.transform.rotate(value=(-math.pi * 0.5), orient_axis='X')  # value = Angle
    bpy.ops.mesh.primitive_plane_add(size=200.0, calc_uvs=True, enter_editmode=False,
                                     align='WORLD', location=(0.0, 0.0, -0.1))

    #======== create a camera and settings
    bpy.ops.object.camera_add(enter_editmode=False, align='WORLD', scale=(5.0, 5.0, 5.0))
    scene.camera = bpy.context.object
    # Set up rotational camera
    cam = bpy.data.objects["Camera"]
    cam.location = (0, -14.5, 5)
    point_at(cam, (0, -14.5, 0.25), roll=math.radians(0))

    scene.cycles.device = 'GPU'

    prefs = bpy.context.preferences
    cprefs = prefs.addons['cycles'].preferences

    # Attempt to set GPU device types if available
    for compute_device_type in ('CUDA', 'OPENCL', 'NONE'):
        try:
            cprefs.compute_device_type = compute_device_type
            break
        except TypeError:
            pass

    # Enable all CPU and GPU devices
    cprefs.get_devices()
    for device in cprefs.devices:
        device.use = True

    #======== create light datablock, set attributes
    light_data = bpy.data.lights.new(name="light_2.80", type='POINT')
    light_data.energy = 12000
    # create new object with our light datablock
    light_object = bpy.data.objects.new(name="light_2.80", object_data=light_data)
    # link light object
    bpy.context.collection.objects.link(light_object)
    # make it active
    bpy.context.view_layer.objects.active = light_object
    # change location
    light_object.location = (-10, -18, 5)

    #======== create another light datablock, set attributes
    light_data1 = bpy.data.lights.new(name="light_top", type='POINT')
    light_data1.energy = 5000
    # # create new object with our light datablock
    light_object1 = bpy.data.objects.new(name="light_top", object_data=light_data1)
    # # link light object
    bpy.context.collection.objects.link(light_object1)
    # # make it active
    bpy.context.view_layer.objects.active = light_object1
    # # change location
    light_object1.location = ( 10, -2, 5 )

    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.device = 'GPU'
    bpy.context.scene.render.resolution_percentage = 100
    bpy.context.scene.cycles.samples = 256
    bpy.context.scene.render.resolution_x = res_x
    bpy.context.scene.render.resolution_y = res_y
    # filename padded with zero
    bpy.context.scene.render.filepath = image_dir + str(i).zfill(4) + ".png"

    # bpy.context.scene.render.filepath = image_dir + str(i).fil + ".png"
    #bpy.context.scene.render.image_settings.compression = 50
    bpy.context.scene.render.image_settings.color_mode = 'RGBA'
    bpy.context.scene.render.image_settings.file_format = 'PNG'
    bpy.ops.render.render(write_still=True)
