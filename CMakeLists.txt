cmake_minimum_required(VERSION 3.25)
project(ITRTG_Generator_Event_Non_Linear_Optimization)

# Set C++ standard globally
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable /Zc:__cplusplus for MSVC to correctly report C++ version
if(MSVC)
    add_compile_options(/Zc:__cplusplus /std:c++20)
endif()

# Optimization flags per compiler
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto -fno-exceptions")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob3 /GL /GS- /Gy /Zi")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    string(REPLACE "/profile" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    string(REPLACE "/PROFILE" "" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif()

# Executables
add_executable(Valentine_2026 "Valentine_2026.cpp")

# Ensure C++23 is applied to targets
set_property(TARGET Valentine_2026 PROPERTY CXX_STANDARD 23)