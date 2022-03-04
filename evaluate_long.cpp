#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <filesystem>
namespace fs = std::filesystem;
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkAbstractArray.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
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

using Vec3f = vtkm::Vec<vtkm::FloatDefault, 3>;

int main(int argc, char **argv)
{
    // load all vtk filenames into a list 
    std::vector<std::string> files;
    std::string path = "/Users/hanmj/Desktop/data/Jet4/";
    for (const auto & entry : fs::directory_iterator(path))
        files.push_back(entry.path());
    std::sort(files.begin(), files.end());

    int start = 0;
    int end = 30;
    int num_seeds = 20;
    int interval = 30;
    int num_fm = (end - start) / interval;
    vtkm::Float32 step_size = 0.01;
    // vec3i seed_dims = vec3i(25, 25, 25);
    vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    double bounds[6];
    reader->SetFileName(files[0].c_str());
    reader->Update();	
    vtkSmartPointer<vtkStructuredPoints> mesh = vtkSmartPointer<vtkStructuredPoints>::New();
    mesh = reader->GetOutput();
    mesh->GetBounds(bounds); 
    std::cout << "Number of points: " << mesh->GetNumberOfPoints() << std::endl;
    vec2f x_range = vec2f(bounds[0], bounds[1]);
    vec2f y_range = vec2f(bounds[2], bounds[3]);
    vec2f z_range = vec2f(bounds[4], bounds[5]);
    std::vector<vec3f> seeds = place_random_seeds_3d_sobol(x_range, y_range, z_range, num_seeds);
    std::cout << "Place " << seeds.size() << std::endl;
    num_seeds = seeds.size();

    vtkm::cont::ArrayHandle<vtkm::Particle> seeds_current;
    seeds_current.Allocate(num_seeds);
    std::vector<std::vector<vec3f>> all_fm;
    all_fm.push_back(seeds);

    for(int f = 0; f < num_fm; f++){
        /* OUTPUT DATA */ 
        std::vector<vtkm::Vec3f_32> pointCoordinates;
        std::vector<vtkm::UInt8> shapes;
        std::vector<vtkm::IdComponent> numIndices;
        std::vector<vtkm::Id> connectivity;
        // Initialize start seeds if f == 0
        if (f == 0){
            for(int i = 0; i < num_seeds; i++){
                vec3f pt = seeds[i];
                seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt.x), static_cast<vtkm::FloatDefault>(pt.y), static_cast<vtkm::FloatDefault>(pt.z)), i));	
                // pointCoordinates.push_back(vtkm::Vec3f_32(pt.x, pt.y, pt.z));
                // shapes.push_back(vtkm::CELL_SHAPE_POLY_LINE);		
                // numIndices.push_back(2);
            }
        }
        // Trace
        for(int c = f * interval; c < (f+1) * interval; c++){
            std::cout << "Starting fm: " << f << std::endl;
            std::cout <<  "Input file " << files[c] << std::endl;
            reader->SetFileName(files[c].c_str());
            reader->Update();	
            vtkSmartPointer<vtkStructuredPoints> mesh = vtkSmartPointer<vtkStructuredPoints>::New();
            mesh = reader->GetOutput();
            int dims[3];
            mesh->GetDimensions(dims);
            int num_pts = mesh->GetNumberOfPoints();
            vtkDataArray* veclocity = mesh->GetPointData()->GetArray(0);
            vtkFloatArray* att1 = vtkFloatArray::SafeDownCast(veclocity);
            vtkm::cont::DataSet dataset;
            vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> velocity_field;
            velocity_field.Allocate(num_pts);

            vtkm::cont::CellSetStructured<3> cellSet;
            cellSet.SetPointDimensions(vtkm::Id3(dims[0], dims[1], dims[2]));
            dataset.SetCellSet(cellSet);
            vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> coords;
            coords.Allocate(num_pts);
            for(int i = 0; i < num_pts; i++){
                double pt[3];
                mesh->GetPoint(i, pt);	
                coords.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])));	
                double vec[3];
                veclocity->GetTuple(i, vec);
                velocity_field.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(vec[0], vec[1], vec[2]));
            }
            dataset.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords));
            dataset.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords));

            const vtkm::cont::DynamicCellSet& cells3d = dataset.GetCellSet();
            const vtkm::cont::CoordinateSystem& coords3d = dataset.GetCoordinateSystem();
    
            // using FieldHandle3d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>>;
            // using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle3d>;
            // using GridEvalType3d = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
            // using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType3d>;
        	using FieldHandle3d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>>;
            using GridEvalType3d = vtkm::worklet::particleadvection::GridEvaluator<FieldHandle3d>;
	        using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType3d>;
	
            // FieldType velocities(velocity_field);
            GridEvalType3d eval_velocity(coords3d, cells3d, velocity_field);
            RK4Type rk4(eval_velocity, static_cast<vtkm::Float32>(step_size));
            vtkm::worklet::ParticleAdvection particleadvection;
            vtkm::worklet::ParticleAdvectionResult res_particles;
            // std::cout << "debug 1" << std::endl;

            
            // advection
            res_particles = particleadvection.Run(rk4, seeds_current, 1);

            auto particles = res_particles.Particles;
            // update seeds_current
            for(int i = 0; i < num_seeds; i++){
                auto next = res_particles.Particles.ReadPortal().Get(i).Pos;
                double pt[3];
                pt[0] = next[0];
                pt[1] = next[1];
                pt[2] = next[2];
                seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])), i));	
            }
        }
        // save current seeds
        std::vector<vec3f> current;
        for(int i = 0; i < num_seeds; i++){
            auto pt = seeds_current.ReadPortal().Get(i).Pos;
            current.push_back(vec3f(pt[0], pt[1], pt[2]));
            // pointCoordinates.push_back(vtkm::Vec3f_32(pt[0], pt[1], pt[2]));
        }
        all_fm.push_back(current);
        // save flow maps 
        // for(int i = 0; i < num_seeds; i++){
        //     for(int n = 0; n < 2; n++){
        //         connectivity.push_back((n*num_seeds)+i);
        //     }
        // }
        // vtkm::cont::DataSetBuilderExplicit dataSetBuilder;
        // vtkm::cont::DataSet outputDS = dataSetBuilder.Create(pointCoordinates, shapes, numIndices, connectivity);
        // std::string outputfile = "fm_" + std::to_string(f) + ".vtk";
        // vtkm::io::VTKDataSetWriter writer(outputfile.c_str());
        // writer.WriteDataSet(outputDS);
    }

    // for(int c = start; c < end; c++){
    //     /* OUTPUT DATA */ 
    //     std::vector<vtkm::Vec3f_32> pointCoordinates;
    //     std::vector<vtkm::UInt8> shapes;
    //     std::vector<vtkm::IdComponent> numIndices;
    //     std::vector<vtkm::Id> connectivity;
    //     std::cout << "Starting cycle: " << c << std::endl;
    //     std::cout <<  "Input file " << files[c] << std::endl;
    //     reader->SetFileName(files[c].c_str());
    //     reader->Update();	
    //     vtkSmartPointer<vtkStructuredPoints> mesh = vtkSmartPointer<vtkStructuredPoints>::New();
    //     mesh = reader->GetOutput();
    //     int dims[3];
    //     mesh->GetDimensions(dims);
    //     std::cout << "Number of points: " << mesh->GetNumberOfPoints() << std::endl;
  	// //     std::cout << "Number of cells: " << mesh->GetNumberOfCells() << std::endl;
  	// //     std::cout << "Array Count : " << mesh->GetPointData()->GetNumberOfArrays() << std::endl;
    // //     std::cout << "Array name: " << mesh->GetPointData()->GetArrayName(0) << std::endl;
  	// //     std::cout << "dimensions: " << dims[0] << " , " << dims[1] << " , " << dims[2] << std::endl;
  	//     std::cout << "bounds: " << bounds[0] << " , " << bounds[1] << " , " << bounds[2] << " , " 
	// 	<< bounds[3] << " , " << bounds[4] << " , " << bounds[5] << std::endl;
        
    //     int num_pts = mesh->GetNumberOfPoints();
    //     vtkDataArray* veclocity = mesh->GetPointData()->GetArray(0);
    //     vtkFloatArray* att1 = vtkFloatArray::SafeDownCast(veclocity);
    //     // if (c == 0){
    //     // Reset seeds_current to start seeds
    //     if ((c+1) % interval == 1){
    //         for(int i = 0; i < num_seeds; i++){
    //             // int index = std::rand() % num_pts;
    //             // double pt[3];
    //             // mesh->GetPoint(index, pt);
    //             vec3f pt = seeds[i];
    //             seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt.x), static_cast<vtkm::FloatDefault>(pt.y), static_cast<vtkm::FloatDefault>(pt.z)), i));	
    //             pointCoordinates.push_back(vtkm::Vec3f_32(pt.x, pt.y, pt.z));
    //             shapes.push_back(vtkm::CELL_SHAPE_POLY_LINE);		
    //             numIndices.push_back(2);
    //             // end-start + 1
    //         }
    //     }

    //     // }
    //     vtkm::cont::DataSet dataset;
    //     vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> velocity_field;
	//     velocity_field.Allocate(num_pts);

	//     vtkm::cont::CellSetStructured<3> cellSet;
    //     cellSet.SetPointDimensions(vtkm::Id3(dims[0], dims[1], dims[2]));
	//     dataset.SetCellSet(cellSet);
    //     vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> coords;
    //     coords.Allocate(num_pts);
    //     for(int i = 0; i < num_pts; i++){
    //         double pt[3];
    //         mesh->GetPoint(i, pt);	
    //         coords.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])));	
    //         double vec[3];
    //         veclocity->GetTuple(i, vec);
    //         velocity_field.WritePortal().Set(i, vtkm::Vec<vtkm::FloatDefault, 3>(vec[0], vec[1], vec[2]));
	//     }
    //     dataset.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords));

	//     const vtkm::cont::DynamicCellSet& cells3d = dataset.GetCellSet();
	//     const vtkm::cont::CoordinateSystem& coords3d = dataset.GetCoordinateSystem();
  
	//     using FieldHandle3d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>>;
    //     using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle3d>;
    //     using GridEvalType3d = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
	//     using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType3d>;
    
    //     FieldType velocities(velocity_field);
	//     GridEvalType3d eval_velocity(coords3d, cells3d, velocities);
	//     RK4Type rk4(eval_velocity, static_cast<vtkm::Float32>(step_size));

	//     vtkm::worklet::ParticleAdvection particleadvection;
	//     vtkm::worklet::ParticleAdvectionResult<vtkm::Particle> res_particles;
    //     std::cout << "debug 1" << std::endl;

    //     // advection
    //     res_particles = particleadvection.Run(rk4, seeds_current, 1);

    //     auto particles = res_particles.Particles;
    //     for(int i = 0; i < num_seeds; i++){
	// 	    auto next = res_particles.Particles.ReadPortal().Get(i).Pos;
    //         double pt[3];
    //         pt[0] = next[0];
    //         pt[1] = next[1];
    //         pt[2] = next[2];
    //         // std::cout << "Particle: " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
    //         // pointCoordinates.push_back(vtkm::Vec3f_32(pt[0], pt[1], pt[2]));
    //         seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])), i));	
    //     }
    //     if ((c+1) % interval == 0) {
    //         for(int i = 0; i < num_seeds; i++){
    //             auto next = res_particles.Particles.ReadPortal().Get(i).Pos;
    //             double pt[3];
    //             pt[0] = next[0];
    //             pt[1] = next[1];
    //             pt[2] = next[2];
    //             // std::cout << "Particle: " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
    //             pointCoordinates.push_back(vtkm::Vec3f_32(pt[0], pt[1], pt[2]));
    //             // seeds_current.WritePortal().Set(i, vtkm::Particle(Vec3f(static_cast<vtkm::FloatDefault>(pt[0]), static_cast<vtkm::FloatDefault>(pt[1]), static_cast<vtkm::FloatDefault>(pt[2])), i));	
    //         }
    //         // save flow maps 
    //         for(int i = 0; i < num_seeds; i++){
	//             for(int n = 0; n < 3; n++){
	// 	            connectivity.push_back((n*num_seeds)+i);
	//             }
    //         }
    //         vtkm::cont::DataSetBuilderExplicit dataSetBuilder;
    //         vtkm::cont::DataSet outputDS = dataSetBuilder.Create(pointCoordinates, shapes, numIndices, connectivity);
    //         std::string outputfile = "fm_" + std::to_string(c+1) + ".vtk";
    //         vtkm::io::VTKDataSetWriter writer(outputfile.c_str());
    //         writer.WriteDataSet(outputDS);
    //     }
    //     std::cout << "debug " << "\n";


    // }
    // for(int i = 0; i < num_seeds; i++){
	//     for(int n = 0; n < (end-start+2); n++){
	// 	connectivity.push_back((n*num_seeds)+i);
	//     }
    // }
    // vtkm::cont::DataSetBuilderExplicit dataSetBuilder;
    // vtkm::cont::DataSet outputDS = dataSetBuilder.Create(pointCoordinates, shapes, numIndices, connectivity);
    // std::string outputfile = "test.vtk";
    // vtkm::io::VTKDataSetWriter writer(outputfile.c_str());
    // writer.WriteDataSet(outputDS);
    return 0;
}