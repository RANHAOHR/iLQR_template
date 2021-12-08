cmake_minimum_required(VERSION 3.19)
project(iLQR)
find_package(PythonLibs REQUIRED)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)


set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/../
        )

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

add_library(matrixOperation
        matrixOperation.cpp
        )

add_executable(iLQR main.cpp)

target_link_libraries(iLQR matrixOperation ${PYTHON_LIBRARIES})


target_include_directories(iLQR
        INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

install(TARGETS iLQR
        EXPORT iLQR_export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(EXPORT iLQR_export
        FILE
        iLQRTargets.cmake
        NAMESPACE
        iLQR::
        DESTINATION
        lib/cmake/iLQR)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/iLQR.h
        DESTINATION
        "include/iLQR")

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/iLQRConfig.cmake
        DESTINATION "lib/cmake/iLQR")

