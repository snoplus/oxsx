set(data_src
    data/DataSetGenerator.cpp
    data/Event.cpp
    data/IO.cpp
    data/LazyOXSXDataSet.cpp
    data/ObsSet.cpp
    data/OXSXDataSet.cpp
    data/ROOTNtuple.cpp
    data/ROOTTree.cpp
    PARENT_SCOPE
)

set(data_headers
    DataSet.h
    DataSetGenerator.h
    Event.h
    IO.h
    LazyOXSXDataSet.h
    ObsSet.h
    ObsSet.hh
    ROOTNtuple.h
    ROOTTree.h
)

file(COPY ${data_headers} DESTINATION ${OXSX_INCLUDE_DIR})