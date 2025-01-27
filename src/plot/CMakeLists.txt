set(plot_src
    plot/ROOTMultiPlot.cpp
    PARENT_SCOPE
)

set(plot_headers
    ROOTMultiPlot.h
)

file(COPY ${plot_headers} DESTINATION ${OXSX_INCLUDE_DIR})