set(cut_src
    cut/BoolCut.cpp
    cut/BoxCut.cpp
    cut/CutCollection.cpp
    cut/CutLog.cpp
    cut/LineCut.cpp
    PARENT_SCOPE
)

set(cut_headers
    BoolCut.h
    BoxCut.h
    Cut.h
    CutCollection.h
    CutLog.h
    LineCut.h
)

file(COPY ${cut_headers} DESTINATION ${OXSX_INCLUDE_DIR})