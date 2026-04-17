#########################################################
# CMake code to build the oxsx library
#########################################################
# Find relevant dependency packages
find_package(ROOT CONFIG REQUIRED Minuit2) # ROOT
find_package(Armadillo REQUIRED) # armadillo (for OXSX)
find_package(HDF5 REQUIRED COMPONENTS CXX C) # HDF5 (for OXSX)
find_package(GSL REQUIRED) # gsl (for OXSX)

# Initialise a bunch of empty lists which we'll be using in a moment...
set(oxsx_src_files "")
set(config_src "")
set(constraint_src "")
set(core_src "")
set(count_src "")
set(cut_src "")
set(data_src "")
set(dist_src "")
set(eventSystematic_src "")
set(fitutil_src "")
set(function_src "")
set(histogram_src "")
set(interval_src "")
set(optimise_src "")
set(plot_src "")
set(rand_src "")
set(systematic_src "")
set(teststat_src "")

# Do some work inside each subdirectory to copy the headers over to a central directory,
# and define the src files present in that sub-drectory.
add_subdirectory(config)
add_subdirectory(constraint)
add_subdirectory(core)
add_subdirectory(count)
add_subdirectory(cut)
add_subdirectory(data)
add_subdirectory(dist)
add_subdirectory(eventSystematic)
add_subdirectory(fitutil)
add_subdirectory(function)
add_subdirectory(histogram)
add_subdirectory(interval)
add_subdirectory(optimise)
add_subdirectory(plot)
add_subdirectory(rand)
add_subdirectory(systematic)
add_subdirectory(teststat)

# Combine all of these source file lists into one!
list(APPEND oxsx_src_files ${config_src} ${constraint_src} ${core_src} ${count_src} ${cut_src} ${data_src} ${dist_src} ${eventSystematic_src} ${fitutil_src} ${function_src} ${histogram_src} ${interval_src} ${optimise_src} ${plot_src} ${rand_src} ${systematic_src} ${teststat_src})

# Create the oxsx library from the source files from all of the sub-directories!
message("oxsx_src_files: ${oxsx_src_files}")
message("OXSX_INCLUDE_DIR: ${OXSX_INCLUDE_DIR}")
add_library(oxsx ${oxsx_src_files} ${OXSX_INCLUDE_DIR})
target_include_directories(oxsx PUBLIC ${OXSX_INCLUDE_DIR})

#---Link external dependencies: ----

# ROOT
target_include_directories(oxsx PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(oxsx PUBLIC ${ROOT_LIBRARIES})
target_compile_options(oxsx PUBLIC "SHELL:${ROOT_DEFINITIONS}")

# HDF5
target_include_directories(oxsx PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(oxsx PUBLIC ${HDF5_CXX_LIBRARIES})

# Armadillo (for OXSX)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(oxsx PUBLIC ${ARMADILLO_LIBRARIES})

# GSL (for OXSX)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(oxsx PUBLIC ${GSL_LIBRARIES})