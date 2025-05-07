cmake_minimum_required(VERSION 3.10)
project(ITRTG Generator Events)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

# Enhanced optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto -fno-exceptions")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob3 /GL /GS- /Gy /Zi")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    # Remove profiling flags which can slow down execution
    string(REPLACE "/profile" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    string(REPLACE "/PROFILE" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif()

# Add executable
add_executable(Easter_v2 "Easter_2025_v2.cpp")