/***********************************************************************************************/
/* A quadratic constraint on the value returned by an arbitrary function of fit parameters.    */
/* For log likelihood and Chi-Square tests this is equivalent to a gaussian contraint on the   */
/* the shape function's returned value.                                                        */
/***********************************************************************************************/
#ifndef __OXSX_SHAPE_CONSTRAINT__
#define __OXSX_SHAPE_CONSTRAINT__
#include <ParameterDict.h>
#include <Exceptions.h>
#include <functional>

/*
 * Firstly, create a typedef for a general function that
 * takes in parameters and returns a value: this is a "shape function".
 * Much like for the Shape Sytematic
 */
typedef std::function<double(const ParameterDict & )> ShapeFunc;

class ShapeConstraint
{
public:
    // Default constructor
    ShapeConstraint() {}
    // Main constructor
    ShapeConstraint(const ShapeFunc &shape_func_, double mean_, double sigma_) : fShapeFunc(shape_func_), fShapeMean(mean_), fShapeSigma(sigma_) {}

    // Evaluate constraint at a particular pair of values
    double Evaluate(const ParameterDict& ) const;

private:
    ShapeFunc fShapeFunc;
    ParameterDict fParamDict;
    double fShapeMean;
    double fShapeSigma;
};
#endif
