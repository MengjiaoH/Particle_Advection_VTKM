#include <iostream>
#include <vector>
#include <chrono>
using namespace std::chrono;
#include "place_seeds.h"
#include "writer.h"

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
#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLImageDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
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
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
// #include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/worklet/particleadvection/RK4Integrator.h>
#include <vtkm/cont/VariantArrayHandle.h>
#include <vtkm/cont/ArrayHandleVirtual.h>
#include <vtkm/worklet/lcs/GridMetaData.h>
#include <vtkm/worklet/lcs/LagrangianStructureHelpers.h>
#include <vtkm/worklet/LagrangianStructures.h>
#include <vtkm/worklet/particleadvection/Stepper.h>


using Vec3f = vtkm::Vec<vtkm::FloatDefault, 3>;
using Vec2f = vtkm::Vec<vtkm::FloatDefault, 2>;

int main(int argc, char** argv)
{
    //Read start location
    std::string start_filename = "/home/mengjiao/Desktop/datasets/flow-maps/heated_cylinder/fm_0.txt";
    std::vector<vec3f> start = read_vec3_from_txt(start_filename, true);
    std::cout << "start location " << start.size() << std::endl;
    // Read end txt 
    std::string end_filename = "/home/mengjiao/Desktop/datasets/flow-maps/heated_cylinder/fm_20.txt";
    std::vector<vec3f> end = read_vec3_from_txt(end_filename, true);
    std::cout << "end location " << end.size() << std::endl;
    
    // Data Dimensions
    vec2i dims = vec2i(30, 90);
    // Read data to get the bounds OR calculate by the range of seeds 
    std::string fname = "/home/mengjiao/Desktop/datasets/boussinesq2d_vti/boussinesq2d.vti";  
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    vtkSmartPointer<vtkImageData> mesh = vtkSmartPointer<vtkImageData>::New();
    // vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    // vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    reader->SetFileName(fname.c_str());
    reader->Update();	
    mesh = reader->GetOutput();
    
    vec2f x_range = vec2f(1000, -1000);
    vec2f y_range = vec2f(1000, -1000);;
    for(int i = 0; i < start.size(); ++i){
        vec3f cur_seed = start[i];
        if(cur_seed.x < x_range.x){
            x_range.x = cur_seed.x;
        }
        if(cur_seed.y < y_range.x){
            y_range.x = cur_seed.y;
        }
        if(cur_seed.x > x_range.y){
            x_range.y = cur_seed.x;
        }
        if(cur_seed.y > y_range.y){
            y_range.y = cur_seed.y;
        }
    }
    std::cout << "x range: " << x_range.x << " " << x_range.y << "\n";
    std::cout << "y range: " << y_range.x << " " << y_range.y << "\n";
    float space_sx = (x_range.y - x_range.x)/(dims.x-1);
    float space_sy = (y_range.y - y_range.x)/(dims.y-1);
    float space_sz = 0.f;
    // Number of steps advanced and the step size 
    int num_steps = 100;
    float step_size = 0.01;

    int seed_pts = start.size();
    vtkm::cont::ArrayHandle<vtkm::Vec2f> lcsInputPoints, lcsOutputPoints;
    lcsInputPoints.Allocate(seed_pts);
    lcsOutputPoints.Allocate(seed_pts);
    for(int i = 0; i < seed_pts; i++){
        lcsInputPoints.WritePortal().Set(i, vtkm::Vec2f(static_cast<vtkm::FloatDefault>(start[i].x),
                                                        static_cast<vtkm::FloatDefault>(start[i].y)));
        lcsOutputPoints.WritePortal().Set(i, vtkm::Vec2f(static_cast<vtkm::FloatDefault>(end[i].x),
                                                         static_cast<vtkm::FloatDefault>(end[i].y)));                                          
    }
    vtkm::cont::DataSetBuilderUniform uniformDatasetBuilder2d;
    vtkm::Id2 outputDims(dims.x, dims.y);
    double x_origin;
    double y_origin;
    double z_origin;
    std::cout << "origin: " << x_origin << " " << y_origin << " " << z_origin << "\n";
    mesh->GetOrigin(x_origin, y_origin, z_origin);
    Vec2f out_origin2d(static_cast<vtkm::FloatDefault>(x_range.x),
                       static_cast<vtkm::FloatDefault>(y_range.x));
    Vec2f out_spacing2d(static_cast<vtkm::FloatDefault>(space_sx),
                        static_cast<vtkm::FloatDefault>(space_sy));
    vtkm::cont::DataSet outputmesh;
    outputmesh = uniformDatasetBuilder2d.Create(outputDims, out_origin2d, out_spacing2d);
    const vtkm::cont::DynamicCellSet& out_cells2d = outputmesh.GetCellSet();
    vtkm::FloatDefault advectionTime = num_steps*step_size;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> ftleField;
    using AnalysisType = vtkm::worklet::LagrangianStructures<2>;
    AnalysisType ftleCalculator(advectionTime, out_cells2d);
    vtkm::worklet::DispatcherMapField<AnalysisType> dispatcher(ftleCalculator);
    dispatcher.Invoke(lcsInputPoints, lcsOutputPoints, ftleField);

    outputmesh.AddPointField("ftle", ftleField);
    // Writing the FTLE
    vtkm::io::VTKDataSetWriter writer("/home/mengjiao/Desktop/datasets/FTLE.vtk");
    writer.WriteDataSet(outputmesh);

    for(int i = 0; i < ftleField.GetNumberOfValues(); i++){
        auto field = ftleField.ReadPortal().Get(i);
    }
    return 0;
}
