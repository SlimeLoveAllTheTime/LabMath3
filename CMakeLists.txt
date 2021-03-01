cmake_minimum_required(VERSION 3.17)
project(LabMath3)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(.idea)

add_executable(LabMath3
        .idea/.gitignore
        .idea/workspace.xml
        Euler.cpp
        Euler.h
        FORSYTHE.h
        main.cpp
        rkf45.cpp
        rkf45.h)
