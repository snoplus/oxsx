#include <Heaviside.h>
#include <Exceptions.h>
#include <iostream>

Heaviside::Heaviside(int nDims_)
{
    fNDims = nDims_;
    fStepPosition.resize(nDims_, 0);
    fSidedness.resize(nDims_, IGNORE);
}

void Heaviside::AddConstraint(unsigned dim_, double pos_, Sidedness side_)
{
    if (dim_ >= GetNDims())
        throw NotFoundError("Heaviside:: Tried to add constraint on non-existent dim!");
    if (side_ != MAXIMUM && side_ != MINIMUM)
        throw ValueError("Heaviside:: Tried to add ambiguous constraint - MAXUMUM or MINIMUM");

    fStepPosition[dim_] = pos_;
    fSidedness[dim_] = side_;
}

double
Heaviside::operator()(const std::vector<double> &vals_) const
{
    if (vals_.size() != GetNDims())
        throw DimensionError("Heaviside::operator()");
    // loop over all the dimensions, return 0 if out of bounds otherwise return 1
    for (size_t i = 0; i < GetNDims(); i++)
    {
        if (fSidedness.at(i) == IGNORE)
            continue;

        double distToStep = vals_.at(i) - fStepPosition.at(i);
        if (fSidedness.at(i) == MINIMUM)
            distToStep *= -1;

        if (distToStep > 0)
            return 0;
    }

    return 1;
}
