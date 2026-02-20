set(constraint_src
    constraint/ConstraintManager.cpp
    constraint/ShapeInterpConstraint.cpp
    constraint/ShapeConstraint.cpp
    PARENT_SCOPE
)

set(constraint_headers
    BivariateQuadraticConstraint.h
    ConstraintManager.h
    QuadraticConstraint.h
    RatioConstraint.h
    ShapeInterpConstraint.h
    ShapeConstraint.h
)

file(COPY ${constraint_headers} DESTINATION ${OXSX_INCLUDE_DIR})