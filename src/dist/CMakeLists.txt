set(dist_src
    dist/AnalyticED.cpp
    dist/BinnedED.cpp
    dist/BinnedEDGenerator.cpp
    dist/CompositeED.cpp
    dist/DistFiller.cpp
    dist/DistTools.cpp
    dist/SpectralFitDist.cpp
    PARENT_SCOPE
)

set(dist_headers
    AnalyticED.h
    BinnedED.h
    BinnedEDGenerator.h
    CompositeED.h
    DistFiller.h
    DistTools.h
    EventDistribution.h
    SpectralFitDist.h
)

file(COPY ${dist_headers} DESTINATION ${OXSX_INCLUDE_DIR})