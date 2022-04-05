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
    std::string data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/tangaroa3d_vtk/U/U_00.00.vtk";
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();  
    double bounds[6];
    reader->SetFileName(data_dir.c_str());
    reader->Update();	
    mesh = reader->GetOutput();
    mesh->GetBounds(bounds); 
    vec2f x_range = vec2f(bounds[0], bounds[1]);
    vec2f y_range = vec2f(bounds[2], bounds[3]);
    vec2f z_range = vec2f(bounds[4], bounds[5]);
    std::cout << "Bounds x : " << bounds[0] << " " << bounds[1] << std::endl;
    std::cout << "Bounds y : " << bounds[2] << " " << bounds[3] << std::endl;
    std::cout << "Bounds z : " << bounds[4] << " " << bounds[5] << std::endl;
    // Place seeds
    std::vector<vec3f> seeds = place_random_seeds_3d(x_range, y_range, z_range, num_seeds);
    std::cout << "Place " << seeds.size() << std::endl;
    for(int i = 0; i < seeds.size(); ++i){
        std::cout << "seed " << i << ": [" << seeds[i].x << ", " << seeds[i].y  << ", " << seeds[i].z <<  "]" << std::endl;
    }
    num_seeds = seeds.size();
    // load vtk 
    // std::string u_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/boussinesq2d_vtk_u/";
    // std::string v_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/boussinesq2d_vtk_v/";
    std::string u_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/tangaroa3d_vtk/U/";
    std::string v_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/tangaroa3d_vtk/V/";
    std::string w_data_dir = "/Users/hanmj/Desktop/projects/Partical_Advection/datasets/tangaroa3d_vtk/W/";
    std::vector<std::filesystem::path> u_filenames;
    std::vector<std::filesystem::path> v_filenames;
    std::vector<std::filesystem::path> w_filenames;
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
    for (const auto& entry : std::filesystem::directory_iterator{w_data_dir}) {
        if (entry.is_regular_file() && (entry.path().extension() == ".vtk")) {
            w_filenames.push_back(entry.path().filename());
        }
    }
    std::sort(w_filenames.begin(), w_filenames.end(), [](const auto& lhs, const auto& rhs) { return lhs.string() < rhs.string();});
    
    vtkm::cont::ArrayHandle<vtkm::Particle> seeds_current;
    seeds_current.Allocate(num_seeds);
    std::vector<std::vector<vec3f>> all_fm;
    all_fm.push_back(seeds);
    // Set start seeds to seeds_current
    for(int i = 0; i < num_seeds; i++){
        vec3f pt = seeds[i];
        seeds_current.WritePortal().Set(i, vtkm::Particle(vtkm::Vec3f(static_cast<vtkm::FloatDefault>(pt.x), static_cast<vtkm::FloatDefault>(pt.y), static_cast<vtkm::FloatDefault>(pt.z)), i));	
    }

    for(int f = 10; f < 11; ++f){
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
        std::cout << "#Points in the data set: " << num_pts << "\n";

        std::string v_filename = v_data_dir + v_filenames[f].string();
        // std::cout << "filename V: " << v_filename << "\n";
        vtkSmartPointer<vtkUnstructuredGridReader> v_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        vtkSmartPointer<vtkUnstructuredGrid> v_mesh = vtkSmartPointer<vtkUnstructuredGrid>::New(); 
        v_reader->SetFileName(v_filename.c_str());
        v_reader->Update();	
        v_mesh = v_reader->GetOutput(); 
        vtkAbstractArray* a2 = v_mesh->GetPointData()->GetArray("V"); 
        vtkFloatArray* att2 = vtkFloatArray::SafeDownCast(a2);

        std::string w_filename = w_data_dir + w_filenames[f].string();
        // std::cout << "filename V: " << v_filename << "\n";
        vtkSmartPointer<vtkUnstructuredGridReader> w_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
        vtkSmartPointer<vtkUnstructuredGrid> w_mesh = vtkSmartPointer<vtkUnstructuredGrid>::New(); 
        w_reader->SetFileName(w_filename.c_str());
        w_reader->Update();	
        w_mesh = w_reader->GetOutput(); 
        vtkAbstractArray* a3 = w_mesh->GetPointData()->GetArray("W"); 
        vtkFloatArray* att3 = vtkFloatArray::SafeDownCast(a3);

        vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> velocity_field;
        velocity_field.Allocate(num_pts);
        for(int i = 0; i < num_pts; i++){ 
            // For Debug
            // if (i < 10)
            //     std::cout << "uv: " << att1->GetTuple1(i) << " " << att2->GetTuple1(i) << "\n";
            velocity_field.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(0.1,  0.1, 0.1)); 
        }
        // Number of Cells in the data set
        int num_cells = v_mesh -> GetNumberOfCells();
        std::cout << "#Cells in the data set: " << num_cells << std::endl;
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
            // std::cout << num << std::endl;
            for(int j = 0; j < num; ++j){
                connectivity.push_back(*(cellPointIds + j));
            }
            shapes.push_back(v_mesh -> GetCellType(i));
            numIndices.push_back(num);
            if ((i > 15670) && (i < 15674)){
                std::cout << "cell " << i << "\n";
                std::cout << "cell id: " << *(cellPointIds + 0) << " " << *(cellPointIds + 1) << " " << *(cellPointIds + 2) << " " << *(cellPointIds + 3) << "\n";
            }
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

        res_particles = particleadvection.Run(rk4, seeds_current, 100);
        auto particles = res_particles.Particles;
        
        // update seeds_current
        for(int i = 0; i < num_seeds; i++){
            auto next = res_particles.Particles.ReadPortal().Get(i).Pos;
            auto num_steps = res_particles.Particles.ReadPortal().Get(i).NumSteps;
            std::cout << "num steps: " << num_steps << "\n";
            double pt[3];
            pt[0] = next[0];
            pt[1] = next[1];
            pt[2] = next[2];
            std::cout << "new pos: " << i << "[" << pt[0] << ", " << pt[1] << ", " << pt[2] << "]" << std::endl;
            seeds_current.WritePortal().Set(i, vtkm::Particle(vtkm::Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])), i));	
        }
        // Save current seeds
        std::vector<vec3f> current;
        for(int i = 0; i < num_seeds; i++){
            auto pt = seeds_current.ReadPortal().Get(i).Pos;
            current.push_back(vec3f(pt[0], pt[1], pt[2]));
        }
        all_fm.push_back(current);
    }

    // write_to_raw(all_fm);
    
    return 0;   
}