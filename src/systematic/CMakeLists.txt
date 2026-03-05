set(systematic_src
    systematic/Convolution.cpp
    systematic/Renorm.cpp
    systematic/Scale.cpp
    systematic/ScaleFunction.cpp
    systematic/Shape.cpp
    systematic/Shift.cpp
    systematic/Systematic.cpp
    PARENT_SCOPE
)

set(systematic_headers
    Convolution.h
    Renorm.h
    Scale.h
    ScaleFunction.h
    Shape.h
    Shift.h
    Systematic.h
)

file(COPY ${systematic_headers} DESTINATION ${OXSX_INCLUDE_DIR})
