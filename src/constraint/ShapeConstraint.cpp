#include "ShapeConstraint.h"

double ShapeConstraint::Evaluate(const ParameterDict& params) const
{
        double val = fShapeFunc(params);
        return (val - fShapeMean) * (val - fShapeMean) / (2 * fShapeSigma * fShapeSigma);
}
