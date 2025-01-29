set(interval_src
    interval/BayesIntervalCalc.cpp
    PARENT_SCOPE
)

set(interval_headers
    BayesIntervalCalc.h
)

file(COPY ${interval_headers} DESTINATION ${OXSX_INCLUDE_DIR})