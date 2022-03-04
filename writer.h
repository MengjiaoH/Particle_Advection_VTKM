#pragma once 
#include <iostream>
#include <fstream>
#include <vector>
#include "place_seeds.h"

void write_to_raw(std::vector<std::vector<vec3f>> fms)
{
    for(int i = 0; i < fms.size(); ++i){
        std::string outfile = "/Users/hanmj/Desktop/projects/Partical_Advection/results/fm_" + std::to_string(i) + ".raw";
        std::ofstream wf(outfile.c_str(), std::ios::out | std::ios::binary);
        if(!wf) {
            std::cout << "Cannot open file!" << std::endl;
        }
        wf.write((char *) &fms[i], sizeof(float) * fms[i].size() * 3);
        wf.close();
        if(!wf.good()) {
            std::cout << "Error occurred at writing time!" << std::endl;
        }
    }
}