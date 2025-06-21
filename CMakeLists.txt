cmake_minimum_required(VERSION 3.16)

project(vpa_stiff_test_2016_arg_arc LANGUAGES CXX)
set(PROJECT_NAME vpa_stiff_test_2016_arg_arc)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core)

add_executable(${PROJECT_NAME}
  main.cpp
  ivp_cauchy_solver.h
  ivp_cauchy_solver.cpp
  csvstream.h
  csvstream.cpp
  newton_iterator.h
  newton_iterator.cpp
  mpTypes.h
  myTimer.h
  config_parser.h
  config_parser.cpp
  debug_mod.h
)
target_link_libraries(${PROJECT_NAME} Qt${QT_VERSION_MAJOR}::Core)

set(EIGEN_ROOT "C:/msys64/mingw64/include/eigen3/")
include_directories(${EIGEN_ROOT})

set(MYCOMPILE_FLAGS "-lmpc -lmpfr -lgmp")
target_compile_options(${PROJECT_NAME} PRIVATE ${MYCOMPILE_FLAGS})

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
