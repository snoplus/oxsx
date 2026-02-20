/***********************************************************************************************/
/* A class to store a set of QuadraticConstraint and BivariateQuadraticConstraint objects,     */
/* associated with fit parameters.                                                             */
/*                                                                                             */
/***********************************************************************************************/
#ifndef __OXSX_CONSTRAINT_MANAGER__
#define __OXSX_CONSTRAINT_MANAGER__
#include <map>
#include <string>
#include <QuadraticConstraint.h>
#include <BivariateQuadraticConstraint.h>
#include <RatioConstraint.h>
#include <FractionalConstraint.h>
#include <ShapeInterpConstraint.h>
#include <Histogram.h>
#include <ParameterDict.h>

class ConstraintManager
{
public:
    // (Default) constructor
    ConstraintManager() : fDebugMode(false) {};
    // Add constraint on a single parameter (symmetric constraint)
    void SetConstraint(const std::string &paramName_, double mean_, double sigma_);
    // Add constraint on a single parameter (asymmetric constraint)
    void SetConstraint(const std::string &paramName_, double mean_, double sigma_hi_, double sigma_lo_);
    // Add correlated constraint between two parameters
    void SetConstraint(const std::string &paramName_1, double mean_1, double sigma_1,
                       const std::string &paramName_2, double mean_2, double sigma_2, double correlation);
    // Add ratio constraint between two parameters
    void SetConstraint(const std::string &paramName_1, const std::string &paramName_2,
                       double ratiomean_, double ratiosigma_);
    // Add shape constraint
    void SetConstraint(const std::string& label, const Histogram& hist);
    // Add fractional diff constraint between two parameters
    void SetFracConstraint(const std::string &paramName_1, const std::string &paramName_2,
                           double fracmean_, double fracsigma_);
    // Evaluate sum of all constraints
    double Evaluate(const ParameterDict &params) const;
    //
    void SetDebugMode(bool debug_mode) { fDebugMode = debug_mode; }

private:
    // Stores quadratic constraints for fit parameters
    std::map<std::string, QuadraticConstraint> fConstraintsInd;
    // Stores correlated constraints between pairs of fit parameters
    std::map<std::pair<std::string, std::string>, BivariateQuadraticConstraint> fConstraintsPair;
    // Stores whether a given fit parameter has an individual (true) or pair (false) constraint
    std::map<std::string, bool> fUseIndConstraints;
    // Stores ratio constraints between pairs of fit parameters
    std::map<std::pair<std::string, std::string>, RatioConstraint> fConstraintsRatio;
    // Stores fractional constraints between pairs of fit parameters
    std::map<std::pair<std::string, std::string>, FractionalConstraint> fConstraintsFrac;
    // Stores any shape constraints, linearly interpolating
    std::map<std::string, ShapeInterpConstraint> fShapeInterpConstraints;
    //
    bool fDebugMode;
};

#endif
