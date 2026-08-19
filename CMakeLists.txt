cmake_minimum_required(VERSION 3.16)
project(pps LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Geant4 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS ROOTNTuple ROOTNTupleUtil)

add_executable(pps pps.c++)
target_include_directories(pps PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(pps PRIVATE ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
target_compile_options(pps PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>)
