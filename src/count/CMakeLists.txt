set(count_src
    count/CountingExperiment.cpp
    count/CountingResult.cpp
    PARENT_SCOPE
)

set(count_headers
    CountingExperiment.h
    CountingResult.h
)

file(COPY ${count_headers} DESTINATION ${OXSX_INCLUDE_DIR})