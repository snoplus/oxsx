set(constraint_src
    constraint/ConstraintManager.cpp
    PARENT_SCOPE
)

set(constraint_headers
    BivariateQuadraticConstraint.h
    ConstraintManager.h
    QuadraticConstraint.h
    RatioConstraint.h
)

file(COPY ${constraint_headers} DESTINATION ${OXSX_INCLUDE_DIR})