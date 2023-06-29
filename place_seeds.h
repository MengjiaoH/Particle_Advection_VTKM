#pragma once 
#include <iostream>
#include <fstream>
#include <vector>
#include "time.h"
#include "sobol.cc"

struct vec3d{
    double x;
    double y;
    double z;
    vec3d(double xx, double yy, double zz){
        x = xx;
        y = yy;
        z = zz;
    }
};

struct vec3f{
    float x;
    float y;
    float z;
    vec3f(float xx, float yy, float zz){
        x = xx;
        y = yy;
        z = zz;
    }
};

struct vec2f{
    float x;
    float y;
    vec2f(float xx, float yy){
        x = xx;
        y = yy;
    }
};

struct vec2i{
    int x;
    int y;
    vec2i(int xx, int yy){
        x = xx;
        y = yy;
    }
};

struct vec3i{
    int x;
    int y;
    int z;
    vec3i(int xx, int yy, int zz){
        x = xx;
        y = yy;
        z = zz;
    }
};

std::vector<vec3f> place_random_seeds_3d(vec2f x_range, vec2f y_range, vec2f z_range, int num_seeds)
{
    std::vector<vec3f> seeds(num_seeds, vec3f(0, 0, 0));
    srand((unsigned) time(NULL));
    for(int i = 0; i < num_seeds; i++){
        float x = (x_range.y - x_range.x) * ((float) rand() / RAND_MAX) + x_range.x;
        float y = (y_range.y - y_range.x) * ((float) rand() / RAND_MAX) + y_range.x;
        float z = (z_range.y - z_range.x) * ((float) rand() / RAND_MAX) + z_range.x;
        seeds[i].x = x;
        seeds[i].y = y;
        seeds[i].z = z;
        // std::cout << x << " " << y << "\n";
    }
    return seeds;
}

std::vector<vec2f> place_random_seeds_2d(vec2f x_range, vec2f y_range, int num_seeds)
{
    std::vector<vec2f> seeds(num_seeds, vec2f(0, 0));
    srand((unsigned) time(NULL));
    for(int i = 0; i < num_seeds; i++){
        float x = (x_range.y - x_range.x) * ((float) rand() / RAND_MAX) + x_range.x;
        float y = (y_range.y - y_range.x) * ((float) rand() / RAND_MAX) + y_range.x;
        seeds[i].x = x;
        seeds[i].y = y;
        // std::cout << x << " " << y << "\n";
    }
    return seeds;
}


std::vector<vec3f> place_sobol_seeds_3d(vec2f x_range, vec2f y_range, vec2f z_range, int num_seeds)
{
    std::vector<vec3f> seeds(num_seeds, vec3f(0, 0, 0));
    // srand((unsigned) time(NULL));
    char* file;
    double **P = sobol_points(num_seeds, 3, file); 
    for(int i = 0; i < num_seeds; i++){
        float x = (x_range.y - x_range.x) * P[i][0] + x_range.x;
        float y = (y_range.y - y_range.x) * P[i][1] + y_range.x;
        float z = (z_range.y - z_range.x) * P[i][2]+ z_range.x;
        seeds[i].x = x;
        seeds[i].y = y;
        seeds[i].z = z;
    }
    return seeds;
}

std::vector<vec2f> place_sobol_seeds_2d(vec2f x_range, vec2f y_range, int num_seeds)
{
    std::vector<vec2f> seeds(num_seeds, vec2f(0, 0));
    // srand((unsigned) time(NULL));
    char* file;
    double **P = sobol_points(num_seeds, 2, file); 
    for(int i = 0; i < num_seeds; i++){
        float x = (x_range.y - x_range.x) * P[i][0] + x_range.x;
        float y = (y_range.y - y_range.x) * P[i][1] + y_range.x;
        seeds[i].x = x;
        seeds[i].y = y;
    }
    return seeds;
}

std::vector<vec3f> load_seeds_from_file(std::string seeds_file, int num_seeds)
{
    std::vector<vec3f> seeds(num_seeds, vec3f(0,0,0));
    std::fstream newfile;
    std::vector<float> temp;
    newfile.open(seeds_file.c_str(), std::ios::in);
    if (newfile.is_open()){   //checking whether the file is open
      std::string tp;
      while(getline(newfile, tp)){ //read data from file object and put it into string.
        std::stringstream ss(tp);
        float value;
        while(ss >> value){
            temp.push_back(value);
        } 
      }
      newfile.close(); //close the file object.
   }
   for(int i = 0 ; i < num_seeds; i++){
       vec3f p = vec3f(temp[3*i], temp[3*i + 1], temp[3*i + 2]);
       seeds[i] = p;
   }
   return seeds;
}

std::vector<vec2f> place_uniform_seeds_2d(vec2f x_range, vec2f y_range, vec2i dims)
{
    
    int num_seeds = dims.x * dims.y;
    std::vector<vec2f> seeds(num_seeds, vec2f(0,0));
    float x_interval = (x_range.y - x_range.x) / (dims.x - 1);
    float y_interval = (y_range.y - y_range.x) / (dims.y - 1);
    for(int j = 0; j < dims.y; j++)
	{
		for(int i = 0; i < dims.x; i++)
		{
			int index = dims.x * j + i;
			float x = i * x_interval + x_range.x;
			float y = j * y_interval + y_range.x;
            vec2f p = vec2f(x, y);
            seeds[index] = p;
		}
	}
    return seeds;
}

std::vector<vec3f> place_uniform_seeds_3d(vec2f x_range, vec2f y_range, vec2f z_range, vec3i dims)
{
    
    int num_seeds = dims.x * dims.y * dims.z;
    std::vector<vec3f> seeds(num_seeds, vec3f(0,0,0));
    float x_interval = (x_range.y - x_range.x) / (dims.x - 1);
    float y_interval = (y_range.y - y_range.x) / (dims.y - 1);
    float z_interval = (z_range.y - z_range.x) / (dims.z - 1);
    for (int k = 0; k < dims.z; k++){

        for(int j = 0; j < dims.y; j++)
	    {
            for(int i = 0; i < dims.x; i++)
            {
                int index = dims.x * dims.y * k + dims.x * j + i;
                float x = i * x_interval + x_range.x;
                float y = j * y_interval + y_range.x;
                float z = k * z_interval + z_range.x;
                // k * z_interval + z_range.x;
                vec3f p = vec3f(x, y, z);
                seeds[index] = p;
            }
	    }
    }
    return seeds;
} 

std::vector<vec3f> place_on_boundary_2d(vec2f x_range, vec2f y_range, vec2i sample_dims)
{
    // int num_seeds = sample_dims.x * 2 + sample_dims.y * 2;
    std::vector<vec3f> seeds;
    float x_interval = (x_range.y - x_range.x) / (sample_dims.x - 1);
    float y_interval = (y_range.y - y_range.x) / (sample_dims.y - 1);
    for(int j = 0; j < sample_dims.y; j++)
	{
			// int index = sample_dims.x * j + i;
			float y = j * y_interval + y_range.x;
            // float y = (y_range.y - y_range.x) * ((float) rand() / RAND_MAX) + y_range.x;
            vec3f p0 = vec3f(x_range.x, y, 1.0);
            vec3f p1 = vec3f(x_range.y, y, 1.0);
            seeds.push_back(p0);
            seeds.push_back(p1);
	}
    for(int j = 0; j < sample_dims.x; j++)
	{
			// int index = sample_dims.x * j + i;
			float x = j * x_interval + x_range.x;
            // float x = (x_range.y - x_range.x) * ((float) rand() / RAND_MAX) + x_range.x;
            vec3f p0 = vec3f(x, y_range.x, 1.0);
            vec3f p1 = vec3f(x, y_range.y, 1.0);
            seeds.push_back(p0);
            seeds.push_back(p1);
	}
	
    

    return seeds;

}
