set(fitutil_src
    fitutil/BinnedEDManager.cpp
    fitutil/BinnedEDShrinker.cpp
    fitutil/EDManager.cpp
    fitutil/EventSystematicManager.cpp
    fitutil/SystematicManager.cpp
    PARENT_SCOPE
)

set(fitutil_headers
    BinnedEDManager.h
    BinnedEDShrinker.h
    EDManager.h
    EventSystematicManager.h
    SystematicManager.h
)

file(COPY ${fitutil_headers} DESTINATION ${OXSX_INCLUDE_DIR})