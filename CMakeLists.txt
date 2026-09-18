cmake_minimum_required(VERSION 3.16)
project(sss LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Geant4 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS ROOTNTuple ROOTNTupleUtil)

add_executable(sss sss.c++)
target_include_directories(sss PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(sss PRIVATE ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
target_compile_options(sss PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>)

add_executable(sss-vis sss-vis.c++)
target_include_directories(sss-vis PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(sss-vis PRIVATE ${ROOT_LIBRARIES})
target_compile_options(sss-vis PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>)
