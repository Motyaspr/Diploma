cmake_minimum_required(VERSION 3.17)
project(Diploma)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -O2")

include_directories(common)
include_directories(viterbi)

add_executable(Diploma
        viterbi/viterbi.cpp
        viterbi/reed_muller.h common/common.h)

add_executable(RMLD RMLD/RMLD.cpp common/common.h viterbi/reed_muller.h)