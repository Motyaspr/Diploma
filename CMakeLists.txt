cmake_minimum_required(VERSION 3.17)
project(Diploma)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -O2")

include_directories(common)
include_directories(viterbi)
include_directories(polar)

add_executable(Diploma
        viterbi/viterbi.cpp
        viterbi/reed_muller.h common/common.h)

add_executable(RMLD RMLD/RMLD.cpp common/common.h viterbi/reed_muller.h)

add_library(polar_libary STATIC common/common.h Polar/polar_encoder.h Polar/TalVardyListDecoder.h
        Polar/TalVardyListDecoder.cpp viterbi/reed_muller.h)

add_executable(polar Polar/encode.cpp)
target_link_libraries(polar polar_libary)