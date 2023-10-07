#pragma once

#include "cgp/cgp.hpp"
#include "../environment.hpp"

struct verre
{
    cgp::mesh_drawable drawable;
    float h;
    float w;

    cgp::vec3 position;
    cgp::rotation_transform orientation;
    cgp::vec3 orientation_vec;

    cgp::vec3 velocity;
    cgp::vec3 force;
    cgp::quaternion angle_velocity;
    cgp::quaternion angle_force;

    void init();
    void draw();
    bool on_floor = false;
};


