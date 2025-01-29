add_subdirectory(data)
add_subdirectory(pdf)
add_subdirectory(plot)

# Define where the example macro files are
set(examples
    beestonBarlowLHH.cpp
    Cuts.cpp
    DataSetGeneration.cpp
    simpleFit.cpp
    simpleSystematicFit.cpp
    systematicFit.cpp
    VaryingCDF.cpp
)

# Create an executable for each macro file
foreach(macro_filename IN LISTS examples)
    message("Macro filename: ${macro_filename}")
    # By stripping macro's filepath, we can obtain string for executable name!
    string(REPLACE ".cpp" "" macro ${macro_filename})
    message("Macro executable name: ${macro}")
    add_executable(${macro} ${macro_filename})
    target_include_directories(${macro} PUBLIC ${OXSX_INCLUDE_DIR})
    target_link_libraries(${macro} PUBLIC oxsx)
endforeach()