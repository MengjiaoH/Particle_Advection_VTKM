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
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/Particles.h>
// #include <vtkm/worklet/particleadvection/RK4Integrator.h>
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
    std::string fname = "/Users/hanmj/Desktop/projects/Partical_Advection/0200.vtk";
    vtkm::Float32 step_size = 0.01;
    int num_seeds = 10;
    // Read Bounds
    vtkSmartPointer<vtkRectilinearGridReader> reader = vtkSmartPointer<vtkRectilinearGridReader>::New();
    double bounds[6];
    reader->SetFileName(fname.c_str());
    reader->Update();	
    vtkSmartPointer<vtkRectilinearGrid> mesh = vtkSmartPointer<vtkRectilinearGrid>::New();
    mesh = reader->GetOutput();
    mesh->GetBounds(bounds); 
    std::cout << "Number of points: " << mesh->GetNumberOfPoints() << std::endl;
    vec2f x_range = vec2f(bounds[0], bounds[1]);
    vec2f y_range = vec2f(bounds[2], bounds[3]);
    vec2f z_range = vec2f(bounds[4], bounds[5]);
    std::cout << "Bounds x : " << bounds[0] << " " << bounds[1] << std::endl;
    std::cout << "Bounds y : " << bounds[2] << " " << bounds[3] << std::endl;
    std::cout << "Bounds z : " << bounds[4] << " " << bounds[5] << std::endl;
    std::vector<vec3f> seeds = place_sobol_seeds_2d(x_range, y_range, num_seeds);
    std::cout << "Place " << seeds.size() << std::endl;
    num_seeds = seeds.size();

    int dims[3] = {512, 512, 1001};
    int interval = 10;
    int num_fm = 100 / 10;
    // (dims[2] - 1) / interval;
    int num_pts = mesh->GetNumberOfPoints();
    std::cout << "There are " << num_pts << " points." << "\n";
    int num_pts_slice = dims[0]*dims[1]; 

    float *xc = (float*) mesh->GetXCoordinates()->GetVoidPointer(0);
    float *yc = (float*) mesh->GetYCoordinates()->GetVoidPointer(0);
    float *zc = (float*) mesh->GetZCoordinates()->GetVoidPointer(0);
    float x_spacing = (xc[dims[0]-1] - xc[0])/(dims[0]-1);
    float y_spacing = (yc[dims[1]-1] - yc[0])/(dims[1]-1);
    float z_spacing = (zc[dims[2]-1] - zc[0])/(dims[2]-1); 
    std::cout << " x spacing: " << x_spacing << "\n";
    std::cout << " y spacing: " << y_spacing << "\n";
    std::cout << " z spacing: " << z_spacing << "\n";
    vtkAbstractArray* a1 = mesh->GetPointData()->GetArray("u");                                   
    vtkAbstractArray* a2 = mesh->GetPointData()->GetArray("v");     
    vtkFloatArray* att1 = vtkFloatArray::SafeDownCast(a1);
    vtkFloatArray* att2 = vtkFloatArray::SafeDownCast(a2);
    
    vtkm::Id3 datasetDims(dims[0], dims[1], dims[2]);
    Vec3f origin3d(static_cast<vtkm::FloatDefault>(xc[0]),
                   static_cast<vtkm::FloatDefault>(yc[0]),
                   static_cast<vtkm::FloatDefault>(zc[0]));
    Vec3f spacing3d(static_cast<vtkm::FloatDefault>(x_spacing),
                    static_cast<vtkm::FloatDefault>(y_spacing),
                    static_cast<vtkm::FloatDefault>(z_spacing));
   
    vtkm::cont::DataSet dataset;
    vtkm::cont::DataSetBuilderUniform uniformDatasetBuilder3d;
    dataset = uniformDatasetBuilder3d.Create(datasetDims, origin3d, spacing3d); 
    
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> velocity_field;
    velocity_field.Allocate(num_pts);
    for(int i = 0; i < num_pts; i++){   
        // std::cout << att1->GetTuple1(i) <<  " " << att2->GetTuple1(i) << "\n";
        velocity_field.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(att1->GetTuple1(i),  att2->GetTuple1(i), 1.0)); 
    }
    

    using FieldHandle3d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>>;     
    const vtkm::cont::DynamicCellSet& cells3d = dataset.GetCellSet(); 
    const vtkm::cont::CoordinateSystem& coords3d = dataset.GetCoordinateSystem(); 
    using GridEvalType3d = vtkm::worklet::particleadvection::GridEvaluator<FieldHandle3d>;
    using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType3d>;
    GridEvalType3d eval_flowmap(coords3d, cells3d, velocity_field);
    RK4Type rk4(eval_flowmap, static_cast<vtkm::Float32>(step_size));
    
    vtkm::worklet::ParticleAdvection particleadvection;
    vtkm::worklet::ParticleAdvectionResult res_particles;

    vtkm::cont::ArrayHandle<vtkm::Particle> seeds_current;
    seeds_current.Allocate(num_seeds);
    std::vector<std::vector<vec3f>> all_fm;
    all_fm.push_back(seeds);

    // Set start seeds to seeds_current
    for(int i = 0; i < num_seeds; i++){
        vec3f pt = seeds[i];
        seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt.x), static_cast<vtkm::FloatDefault>(pt.y), static_cast<vtkm::FloatDefault>(pt.z)), i));	
        // pointCoordinates.push_back(vtkm::Vec3f_32(pt.x, pt.y, pt.z));
        // shapes.push_back(vtkm::CELL_SHAPE_POLY_LINE);		
        // numIndices.push_back(2);
    }
    // std::cout << "debug" << std::endl;
    
    for(int f = 0; f < num_fm; f++){
        // advection
        int num_steps = interval * (f+1);
        res_particles = particleadvection.Run(rk4, seeds_current, interval);
        auto particles = res_particles.Particles;
        // update seeds_current
        for(int i = 0; i < num_pts; i++){
            auto next = res_particles.Particles.ReadPortal().Get(i).Pos;
            double pt[3];
            pt[0] = next[0];
            pt[1] = next[1];
            pt[2] = next[2];
            // std::cout << pt[0] << std::endl;
            seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])), i));	
        }
        // Save current seeds
        std::vector<vec3f> current;
        for(int i = 0; i < num_seeds; i++){
            auto pt = seeds_current.ReadPortal().Get(i).Pos;
            current.push_back(vec3f(pt[0], pt[1], pt[2]));
        }
        all_fm.push_back(current);
    }
    write_to_raw(all_fm);
    return 0;
}