cmake_minimum_required(VERSION 3.17)
project(Diploma)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -O2")

add_executable(Diploma
        viterbi/viterbi.cpp
        CMakeLists.txt viterbi/reed_muller.h)
