set(rand_src
    rand/Rand.cpp
    PARENT_SCOPE
)

set(rand_headers
    Rand.h
)

file(COPY ${rand_headers} DESTINATION ${OXSX_INCLUDE_DIR})