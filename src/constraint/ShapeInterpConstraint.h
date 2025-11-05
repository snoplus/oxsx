/***********************************************************************************************/
/* Sometimes, you have a prior/constraint distribution you want to use within a fit,           */
/* in the form of an n-dimensional grid where each bin has a value corresponding to the        */
/* constraint associated with the bin's centre.                                                */
/* This constraint class performs an n-D linear interpolation algorithm to estimate the        */
/* constraint at any arbitrary point.                                                          */
/***********************************************************************************************/
#ifndef __OXSX_SHAPE_INTERP_CONSTRAINT__
#define __OXSX_SHAPE_INTERP_CONSTRAINT__
#include <Histogram.h>
#include <ParameterDict.h>
#include <vector>


class ShapeInterpConstraint
{
public:
    ShapeInterpConstraint(): fHist() {}
    ShapeInterpConstraint(const Histogram& hist) : fHist(hist) {} // Constructor via Histogram object

    double Evaluate(const ParameterDict& params) const; // Evaluate constraint at a given point in parameter space
private:
    std::vector<double> ExtractConstraintParams(const ParameterDict& params) const;
    double NDLinearInterp(const std::vector<double>& grid_vals, 
                          const std::vector<std::vector<double>>& grid_locs,
                          const std::vector<double>& con_params) const;

    Histogram fHist;
};
#endif