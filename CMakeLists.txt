cmake_minimum_required(VERSION 3.16)
project(Particle_Advection)

# set(CMAKE_C_COMPILER /usr/local/bin/gcc)
# set(CMAKE_CXX_COMPILER /usr/local/bin/g++)
# set(CMAKE_C_FLAGS "-g -fopenmp -Ofast")
set(CMAKE_CXX_STANDARD 17)

# SET(VTKm_DIR "/Users/hanmj/Documents/externals/vtk-m/install/lib/cmake/vtkm-1.5")
SET(VTKm_DIR "/Users/hanmj/Documents/externals/vtk-m-v1.7.1/install/lib/cmake/vtkm-1.7")
find_package(VTKm REQUIRED)

# set(CNPY_DIR "/home/sci/ssane/packages/cnpy/install/include")
# set(CNPY_LIB "/home/sci/ssane/packages/cnpy/install/lib/libcnpy.so")

set(VTK_DIR "/Users/hanmj/Documents/externals/VTK-8.2.0/install/lib/cmake/vtk-8.2")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# add_executable(evaluate_long evaluate_long.cpp)

# target_include_directories(evaluate_long PRIVATE ${VTKm_INCLUDE_DIRS})
#   # ${CNPY_DIR}
# target_link_libraries(evaluate_long PRIVATE vtkm_cont vtkm_io ${VTKm_LIBRARIES} ${VTK_LIBRARIES})
# # ${CNPY_LIB} -fopenmp -lz
# #  target_compile_definitions(ExtractPathlines PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_SERIAL")
# target_compile_options(evaluate_long PRIVATE ${VTKm_COMPILE_OPTIONS})

# add_executable(long_2d long_2d.cpp)

# target_include_directories(long_2d PRIVATE ${VTKm_INCLUDE_DIRS})
# target_link_libraries(long_2d PRIVATE vtkm_cont vtkm_io ${VTKm_LIBRARIES} ${VTK_LIBRARIES})
# # target_compile_definitions(long_2d PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_SERIAL")
# target_compile_options(long_2d PRIVATE ${VTKm_COMPILE_OPTIONS})

add_executable(long_2d_unstructured long_2d_unstructured.cpp)

target_include_directories(long_2d_unstructured PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(long_2d_unstructured PRIVATE vtkm_cont vtkm_io ${VTKm_LIBRARIES} ${VTK_LIBRARIES})
# target_compile_definitions(long_2d PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_SERIAL")
target_compile_options(long_2d_unstructured PRIVATE ${VTKm_COMPILE_OPTIONS})

add_executable(long_3d_unstructured long_3d_unstructured.cpp)

target_include_directories(long_3d_unstructured PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(long_3d_unstructured PRIVATE vtkm_cont vtkm_io ${VTKm_LIBRARIES} ${VTK_LIBRARIES})
# target_compile_definitions(long_2d PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_SERIAL")
target_compile_options(long_3d_unstructured PRIVATE ${VTKm_COMPILE_OPTIONS})