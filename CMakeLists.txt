##

# Author: Yannis STEFANELLI

# Creation Date: 30-05-2023 10:51:18

# Description :

##

cmake_minimum_required(VERSION 3.25.1)

# Executable name
set(PROJECT_NAME
    shader-demo
)

project(${PROJECT_NAME} CXX)

# Source files
file(GLOB SOURCES "src/*.cpp")

if(CMAKE_HOST_SYSTEM MATCHES Linux)
   message(STATUS "Build host runs Linux")

    set( CMAKE_CURRENT_LIB_DIR
        ${CMAKE_CURRENT_SOURCE_DIR}/libs
    )
endif()

if(CMAKE_SYSTEM MATCHES Windows)
   message(STATUS "Target system is Windows")
   set( CMAKE_CURRENT_LIB_DIR
        ${CMAKE_CURRENT_SOURCE_DIR}/windows-libs
   )
endif()

## Libs dir definition
set( SFML_DIR
    ${CMAKE_CURRENT_LIB_DIR}/SFML-2.5.1
)

##

include_directories(
    ${SFML_DIR}/include
)

link_directories(
    ${SFML_DIR}/lib    
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(
    ${PROJECT_NAME}
    sfml-graphics
    sfml-window
    sfml-system
)

add_compile_options(
    -g3
    -Wno-unused-result
)