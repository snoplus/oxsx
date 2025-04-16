set(histogram_src
    histogram/AxisCollection.cpp
    histogram/BinAxis.cpp
    histogram/Histogram.cpp
    histogram/HistTools.cpp
    PARENT_SCOPE
)

set(histogram_headers
    AxisCollection.h
    BinAxis.h
    Histogram.h
    HistTools.h
    HistTools.hpp
)

file(COPY ${histogram_headers} DESTINATION ${OXSX_INCLUDE_DIR})