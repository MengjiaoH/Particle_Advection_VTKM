#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkAbstractArray.h>
#include <vtkRectilinearGridReader.h>
#include <vtkRectilinearGrid.h> 
#include <vtkCellDataToPointData.h>
#include <vtkDataSetWriter.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <stdio.h>
#include <stdlib.h>
#include <vtkCellLocator.h>
#include <vtkm/Matrix.h>
#include <vtkm/Types.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/DynamicCellSet.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ErrorFilterExecution.h>
// #include <vtkm/Particle.h>
// #include <vtkm/filter/flow/worklet/ParticleAdvection.h>
// #include <vtkm/filter/flow/worklet/RK4Integrator.h>
// #include <vtkm/filter/flow/worklet/Stepper.h>
// #include <vtkm/filter/flow/worklet/GridEvaluators.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
// #include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/RK4Integrator.h>
#include <vtkm/worklet/particleadvection/Stepper.h>
#include <vtkm/worklet/particleadvection/Particles.h>

#include <vtkm/cont/VariantArrayHandle.h>
#include <vtkm/cont/ArrayHandleVirtual.h>
#include <vtkm/worklet/lcs/GridMetaData.h>
#include <vtkm/worklet/lcs/LagrangianStructureHelpers.h>
#include <vtkm/worklet/LagrangianStructures.h>


#include "place_seeds.h"
#include "writer.h"

using Vec3f = vtkm::Vec<vtkm::FloatDefault, 3>;

int main(int argc, char **argv)
{
    // load vtk 
    int num_seeds = 500;
    vec2i seed_dims = vec2i(20, 10);
    std::string mode = "uniform";
    // Heated Cylinder
    std::string fname = "/home/mengjiao/Desktop/Examples/Particle_Advection_VTKM/datasets/doublegyre2d_vti/doublegyre2d.vti";
    std::string outfile = "/home/mengjiao/Desktop/Examples/Particle_Advection_VTKM/datasets/flow_maps/double_gyre/short/fm_";
    

    std::string seeds_file = "/home/mengjiao/Desktop/projects/Partical_Advection/result_data/cylinder2d_synthetic/0_1_100/adapted_seeds_45458.txt";
    int interval = 10;
    int start_cycles = 0;
    int stop_cycles = 500;
    
    float offset = 0.0f;
    int num_fm = (stop_cycles - start_cycles) / interval;
    // Heated Cylinder
    vtkm::Float32 step_size = 10.f / 511.f;
    std::cout << "step size: " << step_size << "\n";
    int dims[3] = {256, 128, 512};
    // Cylinder Flow with von Karmen Vortex Street
    // vtkm::Float32 step_size = 0.01;
    // int dims[3] = {640, 80, 1501};
    // Cylinder Flow Synthetic
    // vtkm::Float32 step_size = (5.535 - 1.107) / (500 - 1);
    // int dims[3] = {450, 200, 500};
    // Read Bounds
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    double bounds[6];
    reader->SetFileName(fname.c_str());
    reader->Update();	
    vtkSmartPointer<vtkImageData> mesh = vtkSmartPointer<vtkImageData>::New();
    mesh = reader->GetOutput();
    mesh->GetBounds(bounds); 
    std::cout << "Number of points: " << mesh->GetNumberOfPoints() << std::endl;
    
    vec2f x_range = vec2f(bounds[0] + offset, bounds[1] - offset);
    vec2f y_range = vec2f(bounds[2] + offset, bounds[3] - offset);
    vec2f z_range = vec2f(bounds[4], bounds[5]);
    std::cout << "Bounds x : " << bounds[0] << " " << bounds[1] << std::endl;
    std::cout << "Bounds y : " << bounds[2] << " " << bounds[3] << std::endl;
    std::cout << "Bounds z : " << bounds[4] << " " << bounds[5] << std::endl;
    std::vector<vec3f> seeds;
    if (mode == "uniform"){
        seeds = place_uniform_seeds_2d(x_range, y_range, seed_dims, z_range.x);
    }else if(mode == "random"){
        seeds = place_random_seeds_2d(x_range, y_range, num_seeds, z_range.x);
    }else if(mode == "sobol"){
        seeds = place_sobol_seeds_2d(x_range, y_range, num_seeds, z_range.x);
    }else{
        std::cout << "use seed file" << "\n";
        seeds = load_seeds_from_file(seeds_file, num_seeds, true);
    }
    std::cout << "Place " << seeds.size() << " seeds." << std::endl;
    num_seeds = seeds.size();
    for(int i = 0; i < num_seeds; ++i){
        seeds[i].z = start_cycles * step_size + z_range.x;
    }
    
    int num_pts = mesh->GetNumberOfPoints();
    std::cout << "There are " << num_pts << " points." << "\n";
    int num_pts_slice = dims[0]*dims[1]; 

    double x_spacing;
    double y_spacing;
    double z_spacing;
    mesh->GetSpacing(x_spacing, y_spacing, z_spacing);
    std::cout << " x spacing: " << x_spacing << "\n";
    std::cout << " y spacing: " << y_spacing << "\n";
    std::cout << " z spacing: " << z_spacing << "\n";
    
    double x_origin;
    double y_origin;
    double z_origin;
    
    mesh->GetOrigin(x_origin, y_origin, z_origin);
    std::cout << "origin: " << x_origin << " " << y_origin << " " << z_origin << "\n";
    
    vtkm::Id3 datasetDims(dims[0], dims[1], dims[2]);
    Vec3f origin3d(static_cast<vtkm::FloatDefault>(x_origin),
                   static_cast<vtkm::FloatDefault>(y_origin),
                   static_cast<vtkm::FloatDefault>(z_origin));
    Vec3f spacing3d(static_cast<vtkm::FloatDefault>(x_spacing),
                    static_cast<vtkm::FloatDefault>(y_spacing),
                    static_cast<vtkm::FloatDefault>(z_spacing));
   
    vtkm::cont::DataSet dataset;
    vtkm::cont::DataSetBuilderUniform uniformDatasetBuilder3d;
    dataset = uniformDatasetBuilder3d.Create(datasetDims, origin3d, spacing3d); 

    vtkAbstractArray* a1 = mesh->GetPointData()->GetArray("u");                                   
    vtkAbstractArray* a2 = mesh->GetPointData()->GetArray("v");     
    vtkFloatArray* att1 = vtkFloatArray::SafeDownCast(a1);
    vtkFloatArray* att2 = vtkFloatArray::SafeDownCast(a2);
    
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> velocity_field;
    velocity_field.Allocate(num_pts);
    for(int i = 0; i < num_pts; i++){   
        // std::cout << att1->GetTuple1(i) <<  " " << att2->GetTuple1(i) << "\n";
        velocity_field.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(att1->GetTuple1(i),  att2->GetTuple1(i), 1.0)); 
    }
    using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
    using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
    using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
    using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
    using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;
    FieldType velocities(velocity_field);
    GridEvalType gridEval(dataset.GetCoordinateSystem(), dataset.GetCellSet(), velocities);
    Stepper rk4(gridEval, static_cast<vtkm::Float32>(step_size));
    
    vtkm::worklet::ParticleAdvection particleadvection;
    // vtkm::worklet::ParticleAdvectionResult res_particles;
    vtkm::worklet::ParticleAdvectionResult<vtkm::Particle> res_particles;

    vtkm::cont::ArrayHandle<vtkm::Particle> seeds_current;
    seeds_current.Allocate(num_seeds);
    std::vector<std::vector<vec3f>> all_fm;
    all_fm.push_back(seeds);

    // Set start seeds to seeds_current
    // for(int i = 0; i < num_seeds; i++){
    //     vec3f pt = seeds[i];
    //     seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt.x), static_cast<vtkm::FloatDefault>(pt.y), static_cast<vtkm::FloatDefault>(pt.z)), i));	
    // }

    
    for(int f = start_cycles; f < stop_cycles; f++){
        // advection
        std::cout << f << "\n";
        // Set start seeds to seeds_current
        for(int i = 0; i < num_seeds; i++){
            float t = f * step_size;
            // auto next = particles.ReadPortal().Get(i).Pos;
            // double pt[3];
            // pt[0] = next[0];
            // pt[1] = next[1];
            // pt[2] = next[2];
            // if (i == 50){
            //     // std::cout << pt[0] << " " << pt[1] << " " << pt[2] <<  std::endl;
            //     std::cout << seeds_current.ReadPortal().Get(i).Pos[0] <<  std::endl;
            // }
            vec3f pt = seeds[i];
            seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(t)), i));	
        }
        
        res_particles = particleadvection.Run(rk4, seeds_current, 1);
        auto particles = res_particles.Particles;
        
        
        // Save current seeds
        if((f+1) % interval == 0){
            std::vector<vec3f> current;
            for(int i = 0; i < num_seeds; i++){
                auto pt = particles.ReadPortal().Get(i).Pos;
                current.push_back(vec3f(pt[0], pt[1], pt[2]));
                // pointCoordinates.push_back(vtkm::Vec3f_32(pt[0], pt[1], pt[2]));
            }
            all_fm.push_back(current);
        }
    }
    // for(int i = 0; i < num_seeds; i++){
	//     for(int n = 0; n < (num_fm+1); n++){
    //         // std::cout << i << " " << n << " " << (n * num_seeds) + i << "\n";
	// 	    connectivity.push_back((n*num_seeds)+i);
    //     }
    // }
    // std::string outputfile = "heatedcylinder_test.vtk";
    // vtkm::cont::DataSetBuilderExplicit dataSetBuilder;
    // vtkm::cont::DataSet outputDS = dataSetBuilder.Create(pointCoordinates, shapes, numIndices, connectivity);

    // vtkm::io::VTKDataSetWriter writer(outputfile.c_str());
    // writer.WriteDataSet(outputDS);
    write_to_txt_vec3f(all_fm, outfile);
    // void write_to_txt_vec3f(std::vector<std::vector<vec3f>> fms, std::string outfile)
    return 0;
}
