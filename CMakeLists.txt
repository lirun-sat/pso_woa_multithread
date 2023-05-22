
cmake_minimum_required(VERSION          3.10)
cmake_policy(SET CMP0048                NEW) 

project(pso_woa_planning VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# check if you are in the top-level project or not
if (NOT CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    message(
        STATUS 
        "This project has a top-level one called [${CMAKE_PROJECT_NAME}]"
    )
else()
    message(
        STATUS 
        "This project is a top-level one"
    )
endif()


add_executable(
    ${PROJECT_NAME} 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pso_woa_planning.cpp
)

set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

target_include_directories(${PROJECT_NAME} 
    PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Find pthread library
find_package(Threads REQUIRED)

target_link_libraries(${PROJECT_NAME} 
PRIVATE m
PUBLIC Threads::Threads
)



# gdb调试
SET(CMAKE_BUILD_TYPE "Debug")  
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall") 


                           
                         







