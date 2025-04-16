set(core_src
    core/Combinations.cpp
    core/ComponentManager.cpp
    core/ContainerTools.cpp
    core/DenseMatrix.cpp
    core/ParameterManager.cpp
    core/SparseMatrix.cpp
    PARENT_SCOPE
)

set(core_headers
    Combinations.hpp
    Comparison.hpp
    ComponentManager.h
    ContainerParameter.h
    ContainerParameter.hpp
    ContainerTools.hpp
    DenseMatrix.h
    DoubleParameter.h
    Exceptions.h
    FitComponent.h
    FitParameter.h
    Formatter.hpp
    ParameterDict.h
    ParameterManager.h
    ParameterManager.hpp
    SparseMatrix.h
    TypeTraits.hpp
)

file(COPY ${core_headers} DESTINATION ${OXSX_INCLUDE_DIR})