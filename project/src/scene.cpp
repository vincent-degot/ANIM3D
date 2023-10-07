#include "scene.hpp"


using namespace cgp;

void scene_structure::initialize()
{
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.set_rotation_axis_z();
	camera_control.look_at({ 3.0f, 2.0f, 2.0f }, {0,0,0}, {0,0,1});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());

	obstacle_floor.initialize_data_on_gpu(mesh_primitive_quadrangle({ -1.5f,-1.5f,0 }, { -1.5f,1.5f,0 }, { 1.5f,1.5f,0 }, { 1.5f,-1.5f,0 }));
	obstacle_floor.texture.load_and_initialize_texture_2d_on_gpu(project::path+"assets/wood.jpg");
	obstacle_floor.model.translation = { 0,0,constraint.ground_z };
	obstacle_floor.material.texture_settings.two_sided = true;

	cloth_texture.load_and_initialize_texture_2d_on_gpu(project::path + "assets/cloth.jpg");
	initialize_cloth(gui.N_sample_edge);
    v.init();
}

// Compute a new cloth in its initial position (can be called multiple times)
void scene_structure::initialize_cloth(int N_sample)
{
	cloth.initialize(N_sample);
	cloth_drawable.initialize(N_sample);
	cloth_drawable.drawable.texture = cloth_texture;
	cloth_drawable.drawable.material.texture_settings.two_sided = true;

	constraint.fixed_sample.clear();
	//constraint.add_fixed_position(0, 0, cloth);
	//constraint.add_fixed_position(0, N_sample - 1, cloth);
}

vec3 f;
void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	if (gui.display_frame)
		draw(global_frame, environment);


	// Elements of the scene: Obstacles (floor, sphere), and fixed position
	// ***************************************** //
	
	draw(obstacle_floor, environment);
    /*
	for (auto const& c : constraint.fixed_sample)
	{
		sphere_fixed_position.model.translation = c.second.position;
		draw(sphere_fixed_position, environment);
	}*/

	
	// Simulation of the cloth
	// ***************************************** //
	int const N_step = 1; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	for (int k_step = 0; simulation_running == true && k_step < N_step; ++k_step)
	{
		// Update the forces on each particle
		simulation_compute_force(cloth, parameters);

		// One step of numerical integration
		simulation_numerical_integration(cloth, parameters, parameters.dt);

		// Apply the positional (and velocity) constraints
		simulation_apply_constraints(cloth, constraint);

		// Check if the simulation has not diverged - otherwise stop it
		bool const simulation_diverged = simulation_detect_divergence(cloth);
		if (simulation_diverged) {
			std::cout << "\n *** Simulation has diverged ***" << std::endl;
			std::cout << " > The simulation is stoped" << std::endl;
			simulation_running = false;
		}
	}


	// Cloth display
	// ***************************************** //

	// Prepare to display the updated cloth
	cloth.update_normal();        // compute the new normals
	cloth_drawable.update(cloth); // update the positions on the GPU

	// Display the cloth
	draw(cloth_drawable, environment);
	if (gui.display_wireframe)
		draw_wireframe(cloth_drawable, environment);
		
    vec2 id_cosest = get_closest_(v.position.x, v.position.y);
    vec3 closest_f;
    vec2 coord_closest = cloth.position(id_cosest.x, id_cosest.y);
    bool v_on_tissu = norm(coord_closest - vec2(v.position.x, v.position.y)) < 0.1f;
    if (norm(mouse_coord))
    {
        closest_f = cloth.velocity(id_cosest.x, id_cosest.y);
        vec3 axis_r = normalize(cross(vec3(0,0,1), vec3(mouse_coord, 0)));
        simulation_compute_force_v_angular(v, axis_r, norm(closest_f) / 100, v_on_tissu);
    }
    closest_f = cloth.force(id_cosest.x, id_cosest.y);
    closest_f = cloth.velocity(id_cosest.x, id_cosest.y);
    simulation_compute_force_v(v, closest_f, v_on_tissu);
    simulation_numerical_integration_v(v, parameters, parameters.dt);
    simulation_apply_constraints_v(v, constraint);
    v.drawable.model.translation = (v.position);
    v.drawable.model.rotation = v.orientation;
    draw(v.drawable, environment);

}

void scene_structure::display_gui()
{
	bool reset = false;

	ImGui::Text("Display");
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	ImGui::Checkbox("Texture Cloth", &cloth_drawable.drawable.material.texture_settings.active);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::Text("Simulation parameters");
	ImGui::SliderFloat("Time step", &parameters.dt, 0.0001f, 0.02f, "%.4f", 2.0f);
	ImGui::SliderFloat("Stiffness", &parameters.K, 0.2f, 50.0f, "%.3f", 2.0f);
	ImGui::SliderFloat("Damping", &parameters.mu, 1.0f, 30.0f);
	ImGui::SliderFloat("Mass", &parameters.mass_total, 0.2f, 5.0f, "%.3f", 2.0f);

	ImGui::Spacing(); ImGui::Spacing();

	reset |= ImGui::SliderInt("Cloth samples", &gui.N_sample_edge, 4, 80);

	ImGui::Spacing(); ImGui::Spacing();
	reset |= ImGui::Button("Restart");
	if (reset) {
		initialize_cloth(gui.N_sample_edge);
		simulation_running = true;
	}
}

vec2 scene_structure::get_closest_(float x, float y)
{
    float dist_min = -1;
    vec2 ret;
    int const N = cloth.N_samples();
    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            float dist_x = cloth.position(ku, kv)[0] - x;
            float dist_y = cloth.position(ku, kv)[1] - y;
            float dist_actu = dist_x * dist_x + dist_y * dist_y;
            if (dist_min < 0 || dist_min > dist_actu)
            {
                dist_min = dist_actu;
                ret = vec2(ku, kv);
            }
        }
    }
    return ret;
}

vec2 scene_structure::get_mouse_coord()
{
    vec3 pos = camera_control.camera_model.position();
    vec3 dir = normalize(pos);
    float FOV = camera_projection.field_of_view;
    float ratio = camera_projection.aspect_ratio;
    vec2 xy = inputs.mouse.position.current;

    vec3 cam_r = normalize(cross(camera_control.camera_model.axis_of_rotation, dir));
    vec3 cam_up = normalize(cross(dir, cam_r));

    rotation_transform rot_x = rotation_transform::from_axis_angle(cam_up, -FOV * xy.x);
    rotation_transform rot_y = rotation_transform::from_axis_angle(cam_r, FOV * xy.y / ratio);
    dir = rot_x * rot_y * dir;

    float r = pos[2] / dir[2];
    return vec2(pos.x - (r * dir.x), pos.y - (r * dir.y));
}

void scene_structure::mouse_move_event()
{
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
        if (point_fixed)
        {
            mouse_coord = get_mouse_coord();
            cloth.position(fixed_point_id[0], fixed_point_id[1]) = vec3(mouse_coord, 0);
            cloth.velocity(fixed_point_id[0], fixed_point_id[1]) = vec3(0,0,0);
            cloth.force(fixed_point_id[0], fixed_point_id[1]) = vec3(0,0,0);
        }

}
void scene_structure::mouse_click_event()
{
	camera_control.action_mouse_click(environment.camera_view);
    if (!point_fixed && inputs.keyboard.shift && inputs.mouse.click.left && (inputs.mouse.click.last_action == cgp::last_mouse_cursor_action::click_left))
    {
        vec2 m_coord = get_mouse_coord();
        vec2 closest_id = get_closest_(m_coord.x, m_coord.y);
        fixed_point_id = closest_id;
        point_fixed = true;
    }
    if (point_fixed && !(inputs.keyboard.shift && inputs.mouse.click.left))
    {
        point_fixed = false;
    }
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
    if (point_fixed)
    {
        cloth.velocity(fixed_point_id[0], fixed_point_id[1]) = vec3(0,0,0);
        cloth.force(fixed_point_id[0], fixed_point_id[1]) = vec3(0,0,0);
    }
}

