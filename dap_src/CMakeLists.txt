cmake_minimum_required(VERSION 3.2)
project (dap-g VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 11)
execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions(-DVERSION="${PROJECT_VERSION}" -DUSER="${USER}" -DDATE="${DATE}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(GSL REQUIRED)

add_executable(dap-g 
    src/main.cc 
    src/MLR.h src/MLR.cc 
    src/parser.h src/parser.cc 
    src/controller.h src/controller.cc
)

target_link_libraries(dap-g PUBLIC GSL::gsl GSL::gslcblas)

install(TARGETS dap-g RUNTIME DESTINATION bin)