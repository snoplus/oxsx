set(teststat_src
    teststat/BinnedNLLH.cpp
    teststat/ChiSquare.cpp
    teststat/StatisticSum.cpp
    PARENT_SCOPE
)

set(teststat_headers
    BinnedNLLH.h
    ChiSquare.h
    StatisticSum.h
    TestStatistic.h
)

file(COPY ${teststat_headers} DESTINATION ${OXSX_INCLUDE_DIR})