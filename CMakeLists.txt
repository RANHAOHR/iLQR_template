cmake_minimum_required(VERSION 3.19)
project(iLQR)
find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 14)

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

add_library(matrixOperation
        matrixOperation.cpp
        )

add_executable(ILQR main.cpp)

target_link_libraries(ILQR matrixOperation ${PYTHON_LIBRARIES})