set(config_src
    config/ConfigLoader.cpp
    PARENT_SCOPE
)

set(config_headers
    ConfigLoader.h
    ConfigLoader.hh
    ConfigLoader.hpp
    ini.hpp
)

file(COPY ${config_headers} DESTINATION ${OXSX_INCLUDE_DIR})