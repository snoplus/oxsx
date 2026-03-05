set(constraint_src
    constraint/ConstraintManager.cpp
    constraint/ShapeConstraint.cpp
    constraint/ShapeInterpConstraint.cpp
    PARENT_SCOPE
)

set(constraint_headers
    BivariateQuadraticConstraint.h
    ConstraintManager.h
    FractionalConstraint.h
    QuadraticConstraint.h
    RatioConstraint.h
    ShapeConstraint.h
    ShapeInterpConstraint.h
)

file(COPY ${constraint_headers} DESTINATION ${OXSX_INCLUDE_DIR})