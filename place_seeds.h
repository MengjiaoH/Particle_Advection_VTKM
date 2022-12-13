#pragma once 
#include <iostream>
#include <fstream>
#include <vector>
#include "time.h"
#include "rkcommon/math/vec.h"
#include "sobol.h"
using namespace rkcommon::math;

std::vector<vec3f> place_sobol_seeds_2d(vec2f x_range, vec2f y_range, int num_seeds, float z_start);

std::vector<vec3f> place_sobol_seeds_3d(vec2f x_range, vec2f y_range, vec2f z_range, int num_seeds);

std::vector<vec3f> place_random_seeds_2d(vec2f x_range, vec2f y_range, int num_seeds, float z_start);

std::vector<vec3f> place_random_seeds_3d(vec2f x_range, vec2f y_range, vec2f z_range, int num_seeds);

std::vector<vec3f> load_seeds_from_file(std::string seeds_file, int num_seeds);

std::vector<vec3f> place_uniform_seeds_2d(vec2f x_range, vec2f y_range, vec2i dims, float z_start);

std::vector<vec3f> place_uniform_seeds_3d(vec2f x_range, vec2f y_range, vec2f z_range, vec3i dims);

std::vector<vec3f> place_on_boundary_2d(vec2f x_range, vec2f y_range, vec2i sample_dims);

std::vector<vec3f> load_seeds_from_file(std::string seeds_file, int num_seeds, bool z_zeros);

std::vector<vec2f> place_sobol_seeds_2d_0(vec2f x_range, vec2f y_range, int num_seeds);