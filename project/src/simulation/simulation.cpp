#include "simulation.hpp"

using namespace cgp;

float PI = 3.141592654f;
float deg_to_rad(float deg)
{
    return deg / (180 / PI);
}

float rad_to_deg(float rad)
{
    return rad * (180 / PI);
}

bool eq(float a, float b)
{
    float eps = 0.5f;
    return abs(a - b) < eps;
}

void simulation_compute_force_v_angular(verre& v, vec3 axis_rota, float rota_f, bool v_on_tissu)
{
    const float v_h = v.h;
    const float v_w = v.w;


    quaternion& angle_force = v.angle_force;

    angle_force = rotation_transform::from_axis_angle(vec3(1,0,0), 0.0f).data;


    const vec3 v_ov = v.orientation_vec;
    float sin_theta = abs(v_ov.z) / norm(v_ov);
    float theta = asin(sin_theta);
    if (abs(deg_to_rad(90) - theta) > 0.0001f && abs(v_ov.z) > 0.0001f)
    {
        vec3 axis = normalize(cross(vec3(0,0,1), v_ov));

        float d = std::sqrt(v_h * v_h + v_w * v_w);
        float cos_d_a = v_h / d;
        float diag_angle = acos(cos_d_a);
        float angle_total = diag_angle + theta;

        float rot_f = 0.002f;

        if (cos(angle_total) < 0)
                rot_f = -rot_f;

        if (axis[0] < 0)
            rot_f = -rot_f;
        //if (v_ov.x > 0)
        //    rot_f = -rot_f;
        if (v_ov.y > 0)
            rot_f = -rot_f;
        if (v_ov.z > 0)
            rot_f = -rot_f;
        angle_force *= normalize(rotation_transform::from_axis_angle(axis, -rot_f).data);
    }

    if (rota_f > 0.03f && v_on_tissu && v.on_floor)
        angle_force *= rotation_transform::from_axis_angle(axis_rota, rota_f).data;
}

void simulation_compute_force_v(verre& v, vec3 c_f, bool v_on_tissu)
{
    const float m = 0.05f;
    const vec3 g = { 0,0,-9.81f };
    vec3& force = v.force;
    if (v_on_tissu && v.on_floor)
    {
        force = m * g +  c_f;
    }
    else
        force = m * g;
}
void simulation_numerical_integration_v(verre& v, simulation_parameters const& parameters, float dt)
{
    float const m = 0.05f;

    vec3& vel = v.velocity;
    vec3& p = v.position;
    vec3 const& f = v.force;

    rotation_transform& v_o = v.orientation;
    quaternion& v_angular_vel = v.angle_velocity;
    const quaternion v_angular_f = v.angle_force;
    vec3& v_ov = v.orientation_vec;

    if (v.on_floor)
    {
        vel = vel * 0.3f + dt * f * 80;
        p = p + dt * vel;
    }
    else
    {
        vel = vel + dt * f / m;
        p = p + dt * vel;
    }

    vec3 axis;
    float angle;
    rotation_transform::convert_quaternion_to_axis_angle(v_angular_vel, axis, angle);
    angle *= 0.9f;
    v_angular_vel = rotation_transform::from_axis_angle(axis, angle).data;
    v_angular_vel = normalize(v_angular_vel * v_angular_f);
    auto r = rotation_transform::from_quaternion(normalize(v_angular_vel));
    v_ov = r * v_ov;

    float n = norm(vec3(0,0,1) - v_ov);
    if (norm(vec2(v_ov.x, v_ov.y)) > 0.003f)
    {
        float sin_theta_demi = (n/2) / norm(v_ov);
        float theta = -asin(sin_theta_demi) * 2;
        vec3 ax = normalize(cross(vec3(0,0,1), v_ov));
        //if (ax[0] < 0)
            //theta = -theta;
        auto q = rotation_transform::convert_axis_angle_to_quaternion(ax, theta);
        v_o = rotation_transform::from_quaternion(normalize(q));
    }
    else
        v_o = rotation_transform::from_axis_angle(normalize(v_ov), 0);
}

void simulation_apply_constraints_v(verre& v, constraint_structure const& constraint)
{
    const float v_h = v.h;
    const float v_w = v.w;

    vec3& v_p = v.position;

    const vec3 v_ov = v.orientation_vec;
    float sin_theta = abs(v_ov.z) / norm(v_ov);
    float theta = asin(sin_theta);
    float h;

    if (abs(deg_to_rad(90) - theta) < 0.0001f)
        h = v_h/2;
    else
        h = abs(sin_theta * abs(v_h/2 + 1 / (tan(theta) / (v_w/2))));

    if(v_p[2] < h)
    {
        v_p[2] = h + 0.003f;
        v.on_floor = true;
    }
    else
    {
        v.on_floor = false;
    }
}

vec3 spring(vec3 p1, vec3 p2, float K, float L0)
{
    return K * (norm(p1 - p2) - L0) * ((p2 - p1) / norm(p1 - p2));
}
bool in_grid(int x, int y, int max)
{
    return (x >= 0 && x < max && y >= 0 && y < max);
}

// Fill value of force applied on each particle
// - Gravity
// - Drag
// - Spring force
// - Wind force
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters)
{
    // Direct access to the variables
    //  Note: A grid_2D is a structure you can access using its 2d-local index coordinates as grid_2d(k1,k2)
    //   The index corresponding to grid_2d(k1,k2) is k1 + N1*k2, with N1 the first dimension of the grid.
    //   
    grid_2D<vec3>& force = cloth.force;  // Storage for the forces exerted on each vertex

    grid_2D<vec3> const& position = cloth.position;  // Storage for the positions of the vertices
    grid_2D<vec3> const& velocity = cloth.velocity;  // Storage for the normals of the vertices
    grid_2D<vec3> const& normal = cloth.normal;      // Storage for the velocity of the vertices
    

    size_t const N_total = cloth.position.size();       // total number of vertices
    size_t const N = cloth.N_samples();                 // number of vertices in one dimension of the grid

    // Retrieve simulation parameter
    //  The default value of the simulation parameters are defined in simulation.hpp
    float const K = parameters.K;              // spring stifness
    float const m = parameters.mass_total / N_total; // mass of a particle
    float const mu = parameters.mu;            // damping/friction coefficient
    float const	L0 = 1.0f / (N - 1.0f);        // rest length between two direct neighboring particle
    float const Ldiag = std::sqrt(2 * L0 * L0);


    // Gravity
    const vec3 g = { 0,0,-9.81f };
    for (int ku = 0; ku < N; ++ku)
        for (int kv = 0; kv < N; ++kv)
            force(ku, kv) = m * g;

    // Drag (= friction)
    for (int ku = 0; ku < N; ++ku)
        for (int kv = 0; kv < N; ++kv)
            force(ku, kv) += -mu * m * velocity(ku, kv);


    // TO DO: Add spring forces ...
    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            force(ku, kv) += dot(parameters.wind.direction * parameters.wind.magnitude, normal(kv, kv)) * parameters.wind.direction;
            for (int x = -1; x <= 1; x++)
            {
                for (int y = -1; y <= 1; y++)
                {
                    if (x == 0 && y == 0)
                        continue;
                    if (in_grid(ku + x, kv + y, N))
                    {
                        if (x == 0 || y == 0)
                        {
                            force(ku, kv) += spring(position(ku, kv), position(ku + x, kv + y), K, L0);
                            if (in_grid(ku + x * 2, kv + y * 2, N))
                                force(ku, kv) += spring(position(ku, kv), position(ku + x * 2, kv +  y* 2), K, L0 * 2);
                        }
                        else
                            force(ku, kv) += spring(position(ku, kv), position(ku + x, kv + y), K, Ldiag);
                    }
                }
            }
            // ...
            // force(ku,kv) = ... fill here the force exerted by all the springs attached to the vertex at coordinates (ku,kv).
            // 
            // Notes:
            //   - The vertex positions can be accessed as position(ku,kv)
            //   - The neighbors are at position(ku+1,kv), position(ku-1,kv), position(ku,kv+1), etc. when ku+offset is still in the grid dimension.
            //   - You may want to loop over all the neighbors of a vertex to add each contributing force to this vertex
            //   - To void repetitions and limit the need of debuging, it may be a good idea to define a generic function that computes the spring force between two positions given the parameters K and L0
            //   - If the simulation is a bit too slow, you can speed it up in adapting the parameter N_step in scene.cpp that loops over several simulation step between two displays.
        }
    }

}

void simulation_numerical_integration(cloth_structure& cloth, simulation_parameters const& parameters, float dt)
{
    int const N = cloth.N_samples();
    int const N_total = cloth.position.size();
    float const m = parameters.mass_total/ static_cast<float>(N_total);

    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            vec3& v = cloth.velocity(ku, kv);
            vec3& p = cloth.position(ku, kv);
            vec3 const& f = cloth.force(ku, kv);

            // Standard semi-implicit numerical integration
            v = v + dt * f / m;
            p = p + dt * v;
        }
    }
    
}

void simulation_apply_constraints(cloth_structure& cloth, constraint_structure const& constraint)
{
    // Fixed positions of the cloth
    for (auto const& it : constraint.fixed_sample) {
        position_contraint c = it.second;
        cloth.position(c.ku, c.kv) = c.position; // set the position to the fixed one
    }

    // To do: apply external constraints
    // For all vertex:
    //   If vertex is below floor level ...
    int N = cloth.N_samples();
    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            if (cloth.position(ku,kv).z <= constraint.ground_z)
            {
                cloth.position(ku,kv).z = constraint.ground_z + 0.003f;
            }
        }
    }
    //   If vertex is inside collision sphere ...
}



bool simulation_detect_divergence(cloth_structure const& cloth)
{
    bool simulation_diverged = false;
    const size_t N = cloth.position.size();
    for (size_t k = 0; simulation_diverged == false && k < N; ++k)
    {
        const float f = norm(cloth.force.data.at_unsafe(k));
        const vec3& p = cloth.position.data.at_unsafe(k);

        if (std::isnan(f)) // detect NaN in force
        {
            std::cout << "\n **** NaN detected in forces" << std::endl;
            simulation_diverged = true;
        }

        if (f > 600.0f) // detect strong force magnitude
        {
            std::cout << "\n **** Warning : Strong force magnitude detected " << f << " at vertex " << k << " ****" << std::endl;
            simulation_diverged = true;
        }

        if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)) // detect NaN in position
        {
            std::cout << "\n **** NaN detected in positions" << std::endl;
            simulation_diverged = true;
        }
    }

    return simulation_diverged;
}

