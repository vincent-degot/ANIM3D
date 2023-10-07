#include "verre.hpp"

using namespace cgp;

void verre::init()
{
    drawable.clear();
    h = 0.4;
    w = 0.2;
    position = {0,0,1};
    orientation = rotation_transform();
    orientation_vec = vec3(0,0,1);
    velocity = {0,0,0};
    force = {0,0,0};
    angle_velocity = rotation_transform::from_axis_angle(vec3(0,0,1), 0).data;
    angle_force = rotation_transform::from_axis_angle(vec3(0,0,1), 0).data;
    drawable.initialize_data_on_gpu(mesh_primitive_cylinder(0.1f, {0,0,0.2f}, {0,0,-0.2f}));
}
