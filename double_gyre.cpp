#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <filesystem>
#include <chrono>
using namespace std::chrono;
namespace fs = std::filesystem;
#include "place_seeds.h"
#include <tbb/parallel_for.h>
#include "writer.h"
#include "rkcommon/math/vec.h"

using namespace rkcommon::math;

// #ifndef M_PI
// #define M_PI 3.14159265358979323846264338327950288
// #endif

class doublegyre2d
{
public:
	doublegyre2d(double a = 0.1, double eps = 0.25, double omega = M_PI / 5) : A(a), EPS(eps), OMEGA(omega) {}

	// Samples the vector field of the Double Gyre flow
	vec2f sample(double x, double y, double t)
	{
		return vec2f(
			-A * M_PI*sin(M_PI*(EPS*sin(OMEGA*t)*x*x + (1 - 2 * EPS*sin(OMEGA*t))*x))*cos(M_PI*y),
			A * M_PI*(2 * EPS*sin(OMEGA*t)*x - 2 * EPS*sin(OMEGA*t) + 1)* cos(M_PI*(EPS*sin(OMEGA*t)*x*x + (1 - 2 * EPS*sin(OMEGA*t))*x))*sin(M_PI*y));
	}

private:
	double A; 
	double EPS;
	double OMEGA;
};

int main(int argc, char **argv)
{
    // Define parameters
    int stop_cycle = 300; //trace stop at 300 cycles
    int step_interval = 5; //save end location every 5 cycles
    float step_size = 0.01; // step size for one cycle
    int num_fm = stop_cycle / step_interval; // number of flow maps 

    vec2f x_range = vec2f(0, 2); // domain
    vec2f y_range = vec2f(0, 1);
    int num_seeds = 500;
    vec2i dims = vec2i(256, 128); // resolution

    std::string mode = "random";  // seed placement method 
    /* Place Seeds */
    std::vector<vec3f> seeds;
    if (mode == "uniform"){
        seeds = place_uniform_seeds_2d(x_range, y_range, dims, 0);
    }else if(mode == "random"){
        seeds = place_random_seeds_2d(x_range, y_range, num_seeds, 0);
    }else if(mode == "sobol"){
        seeds = place_sobol_seeds_2d(x_range, y_range, num_seeds, 0);
    }
    
    num_seeds = seeds.size();

    // create double gyre simulator
    doublegyre2d simulations = doublegyre2d();
    
    std::vector<std::vector<vec3f>> fms; // flow maps will only save particle locations along the trajectories
    fms.push_back(seeds); // add start seeds
    
    auto start = high_resolution_clock::now();
    for(int f = 0; f < num_fm; f++){
        std::cout << "#FM: " << f << std::endl;
        tbb::parallel_for( tbb::blocked_range<int>(0, num_seeds), [&](tbb::blocked_range<int> r){
            // for each seed
            for(int p = r.begin(); p < r.end(); ++p){
                vec3f current_point = seeds[p];
                float time = f * step_size * step_interval;
                // advect 5 cycles
                for(int c = 0; c < step_interval; c++){
                    vec2f v = simulations.sample(current_point.x, current_point.y, time);
                    current_point.x = current_point.x + step_size * v.x;
                    current_point.y = current_point.y + step_size * v.y;
                    time += step_size;
                }
                // update current location
                seeds[p].x = current_point.x;
                seeds[p].y = current_point.y;
            }
        });
        // add current location to flow maps
        fms.push_back(seeds);
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time: " << duration.count()<< "s" << std::endl;
    std::cout << "num of fm: " << fms.size() << std::endl;
    std::string outdir = "/home/mengjiao/Desktop/datasets/flow-maps/double_gyre/fm_";
    // std::string outdir = "/Users/hanmj/Desktop/End_2_End_Flow_Vis/datasets/double_gyre/";
    // save_fm_as_raw_2d(fms, outdir);
    write_to_txt_vec3f(fms, outdir);
    
    return 0;
}
