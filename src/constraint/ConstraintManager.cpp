#include <ConstraintManager.h>
#include <QuadraticConstraint.h>
#include <Exceptions.h>
#include <BivariateQuadraticConstraint.h>
#include <ParameterDict.h>
#include <iostream>

void ConstraintManager::SetConstraint(const std::string& paramName_, double mean_, double sigma_) {
    /*
     * Add/update a symmetric quadratic constraint for a single parameter; also do a check
     * that we're not double-counting anywhere.
    */
    if (fUseIndConstraints.count(paramName_) != 0) {
        // Contraint already exists for that parameter:
        // Is it an individual or pair constraint?
        if (fUseIndConstraints[paramName_]) {
            // Is an individual constraint - update!
            fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_);
        } else {
            // There's a pair constraint involving this parameter already - oh no!
            const std::string errmsg = "ConstraintManager::SetConstraint(): cannot set an individual constraint for parameter " + paramName_ + ", as a pair constraint for this parameter already exists!";
            throw LogicError(errmsg);
        }
    } else {
        // Constraint doesn't yet exist for this parameter, so we need to create a new one.
        fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_);
        // Also, note that we now have an individual constraint for this parameter
        fUseIndConstraints[paramName_] = true;
    }
}

void ConstraintManager::SetConstraint(const std::string& paramName_, double mean_, double sigma_hi_, double sigma_lo_) {
    /*
     * Add/update an asymmetric quadratic constraint for a single parameter; also do a check
     * that we're not double-counting anywhere.
    */
    if (fUseIndConstraints.count(paramName_) != 0) {
        // Contraint already exists for that parameter:
        // Is it an individual or pair constraint?
        if (fUseIndConstraints[paramName_]) {
            // Is an individual constraint - update!
            fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_hi_, sigma_lo_);
        } else {
            // There's a pair constraint involving this parameter already - oh no!
            const std::string errmsg = "ConstraintManager::SetConstraint(): cannot set an individual constraint for parameter " + paramName_ + ", as at least one pair constraint for this parameter already exists!";
            throw LogicError(errmsg);
        }
    } else {
        // Constraint doesn't yet exist for this parameter, so we need to create a new one.
        fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_hi_, sigma_lo_);
        // Also, note that we now have an individual constraint for this parameter
        fUseIndConstraints[paramName_] = true;
    }
}

void ConstraintManager::SetConstraint(const std::string& paramName_1, double mean_1, double sigma_1, 
                       const std::string& paramName_2, double mean_2, double sigma_2, double correlation) {
    /*
     * Add/update a bivariate quadratic constraint for a pair of parameters; also do a check to
     * make sure that we don't have constraints for those parameters.
    */
    if (fUseIndConstraints.count(paramName_1) != 0 || fUseIndConstraints.count(paramName_2) != 0) {
        // Constraint already exists for at least one of these params
        // Check if pair constraint already exists
        if (fConstraintsPair.count(std::pair<std::string, std::string>(paramName_1, paramName_2)) != 0) {
            // Pair constraint already exists; update!
            fConstraintsPair[std::pair<std::string, std::string>(paramName_1, paramName_2)] = BivariateQuadraticConstraint(mean_1, sigma_1, mean_2, sigma_2, correlation);
        } else {
            // Oh no, you're trying to add a pair constraint to parameters that already have either
            // individual constraints, or different pair constraints with other parameters!
            const std::string errmsg = "ConstraintManager::SetConstraint(): cannot set a pair constraint for parameters " + paramName_1 + " and " + paramName_2 + ", as at least one other constraint for these parameters already exists!";
            throw LogicError(errmsg);
        }
    } else {
        // Constraints do not yet exists for these params; add a pair constraint.
        fConstraintsPair[std::pair<std::string, std::string>(paramName_1, paramName_2)] = BivariateQuadraticConstraint(mean_1, sigma_1, mean_2, sigma_2, correlation);
        fUseIndConstraints[paramName_1] = false;
        fUseIndConstraints[paramName_2] = false;
    }
}


double ConstraintManager::Evaluate(const ParameterDict& params) const {
    /*
     * Evaluate the value of the sum of all constraints held by this object
    */
    double total = 0;
    // Sum over individual param constraints...
    for (const auto &c_pair : fConstraintsInd) {
        const std::string param_name = c_pair.first;
        const double c = c_pair.second.Evaluate(params.at(param_name));
        total += c;
        if (fDebugMode) { std::cout << "Constraint " << param_name << ": Evaluate() = " << c << std::endl; }
    }
    // ...and then sum over the pair constraints.
    for (const auto &c_pair : fConstraintsPair) {
        const std::pair<std::string, std::string> param_pair = c_pair.first;
        const double c = c_pair.second.Evaluate(params.at(param_pair.first), params.at(param_pair.second));
        total += c;
        if (fDebugMode) { std::cout << "Pair constraint (" << param_pair.first << ", " << param_pair.second << "): Evaluate() = " << c << std::endl; }
    }
    
    return total;
}