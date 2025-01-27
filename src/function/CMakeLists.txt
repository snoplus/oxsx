set(function_src
    function/Gaussian.cpp
    function/GaussianFitter.cpp
    function/Heaviside.cpp
    function/JumpPDF.cpp
    function/SquareRootScale.cpp
    function/VaryingCDF.cpp
    PARENT_SCOPE
)

set(function_headers
    ConditionalPDF.h
    Function.h
    Gaussian.h
    Gaussian.hpp
    GaussianFitter.h
    Heaviside.h
    JumpPDF.h
    PDF.h
    SquareRootScale.h
    VaryingCDF.h
)

file(COPY ${function_headers} DESTINATION ${OXSX_INCLUDE_DIR})