#pragma once 
#include <iostream>
#include <fstream>
#include <vector>
#include "rkcommon/math/vec.h"
using namespace rkcommon::math;

void write_to_raw_2d(std::vector<std::vector<vec2f>> fms)
{
    for(int i = 0; i < fms.size(); ++i){
        std::string outfile = "/Users/hanmj/Desktop/projects/Partical_Advection/results/fm_" + std::to_string(i) + ".raw";
        std::ofstream wf(outfile.c_str(), std::ios::out | std::ios::binary);
        if(!wf) {
            std::cout << "Cannot open file!" << std::endl;
        }
        wf.write((char *) &fms[i], sizeof(float) * fms[i].size() * 2);
        wf.close();
        if(!wf.good()) {
            std::cout << "Error occurred at writing time!" << std::endl;
        }
    }
}
void write_to_raw_3d(std::vector<std::vector<vec3f>> fms)
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

void write_to_txt_2d(std::vector<std::vector<vec2f>> fms, std::string outfile)
{
    for(int i = 0; i < fms.size(); ++i){
        std::string outfilename = outfile + std::to_string(i) + ".txt";
        std::ofstream wf(outfilename.c_str(), std::ios::out);
        if(!wf) {
            std::cout << "Cannot open file!" << std::endl;
        }
        std::vector<vec2f> temp = fms[i];
        // wf.write((char *) &fms[i], sizeof(float) * fms[i].size() * 3);
        for(int j = 0; j < temp.size(); j++){
            wf << temp[j].x << " " << temp[j].y << "\n";
        }
        wf.close();
        if(!wf.good()) {
            std::cout << "Error occurred at writing time!" << std::endl;
        }
    }
}
