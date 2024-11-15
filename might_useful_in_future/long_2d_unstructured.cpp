#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <filesystem>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkAbstractArray.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
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
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayPortalToIterators.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ErrorFilterExecution.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
// #include <vtkm/worklet/particleadvection/IntegratorBase.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/cont/VariantArrayHandle.h>
#include <vtkm/cont/ArrayHandleVirtual.h>
#include <vtkm/worklet/lcs/GridMetaData.h>
#include <vtkm/worklet/lcs/LagrangianStructureHelpers.h>
#include <vtkm/worklet/LagrangianStructures.h>
#include <vtkm/io/VTKDataSetWriter.h>

// Version 1.5
// #include <vtkm/worklet/particleadvection/Integrators.h>

// Version 1.7
#include <vtkm/worklet/particleadvection/RK4Integrator.h>
#include <vtkm/worklet/particleadvection/Stepper.h>
// #include <vtkm/worklet/testing/GenerateTestDataSets.h>


#include "place_seeds.h"
#include "writer.h"

using Vec2f = vtkm::Vec<vtkm::FloatDefault, 2>;

int main(int argc, char **argv)
{
    vtkm::Float32 step_size = 0.01;
    int num_seeds = 10;
    std::string mode = "uniform";
    vec2i seed_dims = vec2i(5, 10);
    std::string outfile = "/Users/hanmj/Desktop/projects/Partical_Advection/results/fm_";
    std::string data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/boussinesq2d_vtk_u/U_00.00.vtk";
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();  
    double bounds[4];
    reader->SetFileName(data_dir.c_str());
    reader->Update();	
    mesh = reader->GetOutput();
    mesh->GetBounds(bounds); 
    vec2f x_range = vec2f(bounds[0], bounds[1]);
    vec2f y_range = vec2f(bounds[2], bounds[3]);
    std::cout << "Bounds x : " << bounds[0] << " " << bounds[1] << std::endl;
    std::cout << "Bounds y : " << bounds[2] << " " << bounds[3] << std::endl;
    // Place seeds
    std::vector<vec2f> seeds;
    if (mode == "uniform"){
        seeds = place_uniform_seeds_2d(x_range, y_range, seed_dims);
    }else if(mode == "sobol"){
        seeds = place_sobol_seeds_2d(x_range, y_range, num_seeds);
    }else if(mode == "random"){
        seeds = place_random_seeds_2d(x_range, y_range, num_seeds);
    }
    
    std::cout << "Place " << seeds.size() << std::endl;
    for(int i = 0; i < seeds.size(); ++i){
        std::cout << "seed " << i << ": [" << seeds[i].x << ", " << seeds[i].y  << "]" << std::endl;
    }
    num_seeds = seeds.size();
    // load vtk 
    // std::string u_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/boussinesq2d_vtk_u/";
    // std::string v_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/boussinesq2d_vtk_v/";
    std::string u_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/subset_u/";
    std::string v_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/subset_v/";
    std::vector<std::filesystem::path> u_filenames;
    std::vector<std::filesystem::path> v_filenames;
    for (const auto& entry : std::filesystem::directory_iterator{u_data_dir}) {
        if (entry.is_regular_file() && (entry.path().extension() == ".vtk")) {
            u_filenames.push_back(entry.path().filename());
        }
    }
    std::sort(u_filenames.begin(), u_filenames.end(), [](const auto& lhs, const auto& rhs) { return lhs.string() < rhs.string();});
        for (const auto& entry : std::filesystem::directory_iterator{v_data_dir}) {
        if (entry.is_regular_file() && (entry.path().extension() == ".vtk")) {
            v_filenames.push_back(entry.path().filename());
        }
    }
    std::sort(v_filenames.begin(), v_filenames.end(), [](const auto& lhs, const auto& rhs) { return lhs.string() < rhs.string();});
    vtkm::cont::ArrayHandle<vtkm::Particle> seeds_current;
    seeds_current.Allocate(num_seeds);
    std::vector<std::vector<vec2f>> all_fm;
    all_fm.push_back(seeds);
    // Set start seeds to seeds_current
    for(int i = 0; i < num_seeds; i++){
        vec2f pt = seeds[i];
        seeds_current.WritePortal().Set(i, vtkm::Particle(vtkm::Vec3f(static_cast<vtkm::FloatDefault>(pt.x), static_cast<vtkm::FloatDefault>(pt.y), 0), i));	
    }

    for(int f = 0; f < u_filenames.size(); ++f){
        /* ************************************************************* */
        // Load velocity u and v 
        /* ************************************************************* */
        std::string u_filename = u_data_dir + u_filenames[f].string();
        // std::cout << "filename U: " << u_filename << "\n";
        vtkSmartPointer<vtkUnstructuredGridReader> u_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        vtkSmartPointer<vtkUnstructuredGrid> u_mesh = vtkSmartPointer<vtkUnstructuredGrid>::New(); 
        u_reader->SetFileName(u_filename.c_str());
        u_reader->Update();	
        u_mesh = u_reader->GetOutput(); 
        vtkAbstractArray* a1 = u_mesh->GetPointData()->GetArray("U"); 
        vtkFloatArray* att1 = vtkFloatArray::SafeDownCast(a1);
        // Number of Points in the  data set 
        int num_pts = u_mesh->GetNumberOfPoints();
        // std::cout << "#Points in the data set: " << num_pts << "\n";

        std::string v_filename = v_data_dir + v_filenames[f].string();
        std::cout << "filename V: " << v_filename << "\n";
        vtkSmartPointer<vtkUnstructuredGridReader> v_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        vtkSmartPointer<vtkUnstructuredGrid> v_mesh = vtkSmartPointer<vtkUnstructuredGrid>::New(); 
        v_reader->SetFileName(v_filename.c_str());
        v_reader->Update();	
        v_mesh = v_reader->GetOutput(); 
        vtkAbstractArray* a2 = v_mesh->GetPointData()->GetArray("V"); 
        vtkFloatArray* att2 = vtkFloatArray::SafeDownCast(a2);
        vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> velocity_field;
        velocity_field.Allocate(num_pts);
        for(int i = 0; i < num_pts; i++){ 
            // For Debug
            // if (i < 10)
            //     std::cout << "uv: " << att1->GetTuple1(i) << " " << att2->GetTuple1(i) << "\n";
            velocity_field.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(att1->GetTuple1(i),  att2->GetTuple1(i), 0)); 
        }
        // Number of Cells in the data set
        int num_cells = v_mesh -> GetNumberOfCells();
        // std::cout << "#Cells in the data set: " << num_cells << std::endl;
        /* ************************************************************* */
        // Create dataset 
        /* ************************************************************* */
        std::vector<vtkm::Vec3f_32> pointCoordinates;
        std::vector<vtkm::Float32> pointCoord_x;
        std::vector<vtkm::Float32> pointCoord_y;
        std::vector<vtkm::UInt8> shapes;
        std::vector<vtkm::IdComponent> numIndices;
        std::vector<vtkm::Id> connectivity;
        for(int i = 0; i < num_cells; ++i){
            vtkIdType num = 0;
            vtkIdType* cellPointIds;
            v_mesh ->GetCellPoints(i, num, cellPointIds);
            connectivity.push_back(*(cellPointIds + 0));
            connectivity.push_back(*(cellPointIds + 1));
            connectivity.push_back(*(cellPointIds + 3));
            connectivity.push_back(*(cellPointIds + 2));
            // for(int j = 0; j < num; ++j){
            //     connectivity.push_back(*(cellPointIds + j));
            // }
            shapes.push_back(vtkm::UInt8(9));
            // v_mesh -> GetCellType(i)
            numIndices.push_back(num);
        }
        
        for(int i = 0; i < num_pts; ++i){
            double* point = v_mesh ->GetPoint(i);
            // std::cout << point[0] << " " << point[1] <<  " " << point[2] << "\n";
            pointCoordinates.push_back(vtkm::Vec3f(point[0], point[1], point[2]));
            // pointCoord_x.push_back(*(point+0));
            // pointCoord_y.push_back(*(point+1));
        }

        vtkm::cont::DataSetBuilderExplicit dataSetBuilder;
        vtkm::cont::DataSet dataset = dataSetBuilder.Create(pointCoordinates, shapes, numIndices, connectivity); // 3D 
        // vtkm::cont::DataSet dataset = dataSetBuilder.Create(pointCoord_x, pointCoord_y, shapes, numIndices, connectivity); // 2D
        vtkm::io::VTKDataSetWriter writer("test.vtk");
        writer.WriteDataSet(dataset);
        // std::cout << "Debug boundx: " << bound.X << std::endl;
        // std::cout << "Debug boundy: " << bound.Y << std::endl;
        // // std::cout << "Debug boundz: " << bound.Z << std::endl;
        // std::cout << "Debug num cells: " << cells3d.GetNumberOfCells() << "\n";
        // Version 1.5
        // using FieldHandle3d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>>;     
        // const vtkm::cont::DynamicCellSet& cells3d = dataset.GetCellSet(); 
        // const vtkm::cont::CoordinateSystem& coords3d = dataset.GetCoordinateSystem(); 
        // vtkm::Bounds bound = coords3d.GetBounds();
        // using GridEvalType3d = vtkm::worklet::particleadvection::GridEvaluator<FieldHandle3d>;
        // using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType3d>;
        // GridEvalType3d eval_flowmap(coords3d, cells3d, velocity_field);
        // RK4Type rk4(eval_flowmap, static_cast<vtkm::Float32>(step_size));
        // Version 1.7
        using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec3f>;
        using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
        using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
        using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
        using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;
        FieldType velocities(velocity_field);
        GridEvalType gridEval(dataset.GetCoordinateSystem(), dataset.GetCellSet(), velocities);
        Stepper rk4(gridEval, static_cast<vtkm::Float32>(step_size));
        

        vtkm::worklet::ParticleAdvection particleadvection;
        // Version 1.5
        // vtkm::worklet::ParticleAdvectionResult res_particles;
        // Version 1.7
        vtkm::worklet::ParticleAdvectionResult<vtkm::Particle> res_particles;

        res_particles = particleadvection.Run(rk4, seeds_current, 1);
        auto particles = res_particles.Particles;
        
        // update seeds_current
        for(int i = 0; i < num_seeds; i++){
            auto next = res_particles.Particles.ReadPortal().Get(i).Pos;
            // auto num_steps = res_particles.Particles.ReadPortal().Get(i).NumSteps;
            // std::cout << "num steps: " << num_steps << "\n";
            double pt[3];
            pt[0] = next[0];
            pt[1] = next[1];
            pt[2] = next[2];
            // std::cout << "new pos: " << i << "[" << pt[0] << ", " << pt[1] << ", " << pt[2] << "]" << std::endl;
            seeds_current.WritePortal().Set(i, vtkm::Particle(vtkm::Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])), i));	
        }
        // Save current seeds
        std::vector<vec2f> current;
        for(int i = 0; i < num_seeds; i++){
            auto pt = seeds_current.ReadPortal().Get(i).Pos;
            current.push_back(vec2f(pt[0], pt[1]));
        }
        all_fm.push_back(current);
    }
    
    write_to_txt_2d(all_fm, outfile);
    
    return 0;   
}