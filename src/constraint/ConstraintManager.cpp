#include <ConstraintManager.h>

void ConstraintManager::SetConstraint(const std::string &paramName_, double mean_, double sigma_)
{
    /*
     * Add/update a symmetric quadratic constraint for a single parameter; also do a check
     * that we're not double-counting anywhere.
     */
    if (fUseIndConstraints.count(paramName_) != 0)
    {
        // Contraint already exists for that parameter:
        // Is it an individual or pair constraint?
        if (fUseIndConstraints[paramName_])
        {
            // Is an individual constraint - update!
            fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_);
        }
        else
        {
            // There's a pair constraint involving this parameter already - oh no!
            const std::string errmsg = "ConstraintManager::SetConstraint(): cannot set an individual constraint for parameter " + paramName_ + ", as a pair constraint for this parameter already exists!";
            throw LogicError(errmsg);
        }
    }
    else
    {
        // Constraint doesn't yet exist for this parameter, so we need to create a new one.
        fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_);
        // Also, note that we now have an individual constraint for this parameter
        fUseIndConstraints[paramName_] = true;
    }
}

void ConstraintManager::SetConstraint(const std::string &paramName_, double mean_, double sigma_hi_, double sigma_lo_)
{
    /*
     * Add/update an asymmetric quadratic constraint for a single parameter; also do a check
     * that we're not double-counting anywhere.
     */
    if (fUseIndConstraints.count(paramName_) != 0)
    {
        // Contraint already exists for that parameter:
        // Is it an individual or pair constraint?
        if (fUseIndConstraints[paramName_])
        {
            // Is an individual constraint - update!
            fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_hi_, sigma_lo_);
        }
        else
        {
            // There's a pair constraint involving this parameter already - oh no!
            const std::string errmsg = "ConstraintManager::SetConstraint(): cannot set an individual constraint for parameter " + paramName_ + ", as at least one pair constraint for this parameter already exists!";
            throw LogicError(errmsg);
        }
    }
    else
    {
        // Constraint doesn't yet exist for this parameter, so we need to create a new one.
        fConstraintsInd[paramName_] = QuadraticConstraint(mean_, sigma_hi_, sigma_lo_);
        // Also, note that we now have an individual constraint for this parameter
        fUseIndConstraints[paramName_] = true;
    }
}

void ConstraintManager::SetConstraint(const std::string &paramName_1, double mean_1, double sigma_1,
                                      const std::string &paramName_2, double mean_2, double sigma_2, double correlation)
{
    /*
     * Add/update a bivariate quadratic constraint for a pair of parameters; also do a check to
     * make sure that we don't have constraints for those parameters.
     */
    if (fUseIndConstraints.count(paramName_1) != 0 || fUseIndConstraints.count(paramName_2) != 0)
    {
        // Constraint already exists for at least one of these params
        // Check if pair constraint already exists
        if (fConstraintsPair.count(std::pair<std::string, std::string>(paramName_1, paramName_2)) != 0)
        {
            // Pair constraint already exists; update!
            fConstraintsPair[std::pair<std::string, std::string>(paramName_1, paramName_2)] = BivariateQuadraticConstraint(mean_1, sigma_1, mean_2, sigma_2, correlation);
        }
        else
        {
            // Oh no, you're trying to add a pair constraint to parameters that already have either
            // individual constraints, or different pair constraints with other parameters!
            const std::string errmsg = "ConstraintManager::SetConstraint(): cannot set a pair constraint for parameters " + paramName_1 + " and " + paramName_2 + ", as at least one other constraint for these parameters already exists!";
            throw LogicError(errmsg);
        }
    }
    else
    {
        // Constraints do not yet exists for these params; add a pair constraint.
        fConstraintsPair[std::pair<std::string, std::string>(paramName_1, paramName_2)] = BivariateQuadraticConstraint(mean_1, sigma_1, mean_2, sigma_2, correlation);
        fUseIndConstraints[paramName_1] = false;
        fUseIndConstraints[paramName_2] = false;
    }
}

void ConstraintManager::SetConstraint(const std::string &paramName_1, const std::string &paramName_2, double ratiomean_, double ratiosigma_)
{
    /*
     * Add/update a ratio constraint for a pair of parameters; we don't need to worry about other constraints
     */
    fConstraintsRatio[std::pair<std::string, std::string>(paramName_1, paramName_2)] = RatioConstraint(ratiomean_, ratiosigma_);
}

void ConstraintManager::SetConstraint(const ParameterDict &params_, ShapeFunc func_, double mean_, double sigma_)
{
    // Add an analytical shape function
    fShapeConstraints[params_] = ShapeConstraint(func_, mean_, sigma_);
}


void ConstraintManager::SetConstraint(const std::string& label, const Histogram& hist)
{
    /*
     * Add/update an interpolated shape constraint for a number of parameters.
     * Histogram gives the lattice of prior/constraint values at the bins' centres;
     * label is just a string to refer to this constraint, allowing for it to be updated later.
     */
    fShapeInterpConstraints[label] = ShapeInterpConstraint(hist);
}

void ConstraintManager::SetFracConstraint(const std::string &paramName_1, const std::string &paramName_2, double fracmean_, double fracsigma_)
{
    /*
     * Add/update a fractional difference constraint for a pair of parameters; we don't need to worry about other constraints
     */
    fConstraintsFrac[std::pair<std::string, std::string>(paramName_1, paramName_2)] = FractionalConstraint(fracmean_, fracsigma_);
}

double ConstraintManager::Evaluate(const ParameterDict &params) const
{
    /*
     * Evaluate the value of the sum of all constraints held by this object
     */
    double total = 0;
    // Sum over individual param constraints...
    for (const auto &c_pair : fConstraintsInd)
    {
        const std::string param_name = c_pair.first;
        const double c = c_pair.second.Evaluate(params.at(param_name));
        total += c;
        if (fDebugMode)
        {
            std::cout << "Constraint " << param_name << ": Evaluate() = " << c << std::endl;
        }
    }

    // ...and then sum over the pair constraints...
    for (const auto &c_pair : fConstraintsPair)
    {
        const std::pair<std::string, std::string> param_pair = c_pair.first;
        const double c = c_pair.second.Evaluate(params.at(param_pair.first), params.at(param_pair.second));
        total += c;
        if (fDebugMode)
        {
            std::cout << "Pair constraint (" << param_pair.first << ", " << param_pair.second << "): Evaluate() = " << c << std::endl;
        }
    }

    // ...and then sum over the ratio constraints...
    for (const auto &c_pair : fConstraintsRatio)
    {
        const std::pair<std::string, std::string> param_pair = c_pair.first;
        const double c = c_pair.second.Evaluate(params.at(param_pair.first), params.at(param_pair.second));
        total += c;
        if (fDebugMode)
        {
            std::cout << "Ratio constraint (" << param_pair.first << ", " << param_pair.second << "): Evaluate() = " << c << std::endl;
        }
    }

    // ...and then sum over the fractional constraints...
    for (const auto &c_pair : fConstraintsFrac)
    {
        const std::pair<std::string, std::string> param_pair = c_pair.first;
        const double c = c_pair.second.Evaluate(params.at(param_pair.first), params.at(param_pair.second));
        total += c;
        if (fDebugMode)
        {
            std::cout << "Fractional constraint (" << param_pair.first << ", " << param_pair.second << "): Evaluate() = " << c << std::endl;
        }
    }

    // ...and sum over the shape interpolation constraints
    for (const auto &c_pair : fShapeInterpConstraints)
    {
        const double c = c_pair.second.Evaluate(params);
        total += c;
        if (fDebugMode)
        {
            std::cout << "Shape Interp. Constraint " << c_pair.first << ": Evaluate() = " << c << std::endl;
        }
    }

    // ...and finally sum over the analytical shape constraints
    for (const auto &c_pair : fShapeConstraints)
    {
        const double c = c_pair.second.Evaluate(params);
        total += c;
        if (fDebugMode)
        {
            std::cout << "Shape Constraint : Evaluate() = " << c << std::endl;
        }
    }

    return total;
}
