set(eventSystematic_src
    eventSystematic/EventConvolution.cpp
    eventSystematic/EventReconvolution.cpp
    eventSystematic/EventScale.cpp
    eventSystematic/EventShift.cpp
    eventSystematic/EventSystematic.cpp
    PARENT_SCOPE
)

set(eventSystematic_headers
    EventConvolution.h
    EventReconvolution.h
    EventScale.h
    EventShift.h
    EventSystematic.h
)

file(COPY ${eventSystematic_headers} DESTINATION ${OXSX_INCLUDE_DIR})