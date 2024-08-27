#include <iostream>
template <typename T>
void LeapFrog::Hamiltonian(ParameterDict &qs_, ParameterDict &ps_,
                           const ParameterDict &masses_,
                           double epsilon_, T &differ_)
{
    // step 1
    ParameterDict grad = differ_(qs_);
    for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_ / 2 * grad[it->first];

    // step 2
    for (ParameterDict::iterator it = qs_.begin(); it != qs_.end(); ++it)
        it->second += epsilon_ * ps_[it->first] / masses_.at(it->first);

    // step 3
    grad = differ_(qs_);
    for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_ / 2 * grad[it->first];
}

template <typename T>
void LeapFrog::Hamiltonian(ParameterDict &qs_, ParameterDict &ps_,
                           const ParameterDict &masses_,
                           double epsilon_, T &differ_, int nSteps_)
{
    // step 1 take a half step in momentum
    ParameterDict grad = differ_(qs_);
    for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_ / 2 * grad[it->first];

    // step 2 alternate full steps for position and momentum
    for (int i = 0; i < nSteps_; i++)
    {
        // Full step for position
        for (ParameterDict::iterator it = qs_.begin(); it != qs_.end(); ++it)
            it->second += epsilon_ * ps_[it->first] / masses_.at(it->first);

        // Full step for momentum, unless this is the last step
        if (i != nSteps_ - 1)
        {
            grad = differ_(qs_);
            for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
                it->second -= epsilon_ * grad[it->first];
        }
    }

    // step 3 make a half step for momentum
    grad = differ_(qs_);
    for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
        it->second -= epsilon_ / 2 * grad[it->first];
}

template <typename T>
void LeapFrog::Hamiltonian(ParameterDict &qs_, ParameterDict &ps_,
                           const ParameterDict &masses_,
                           double epsilon_, T &differ_, int nSteps_,
                           const ParameterDict &minima_,
                           const ParameterDict &maxima_)
{
    /// https://papers.nips.cc/paper/5801-reflection-refraction-and-hamiltonian-monte-carlo.pdf
    /// pg4 - simplified for only 1d constant boundaries and no refraction (= hard boundaries with delta U = inf)

    for (int i = 0; i < nSteps_; i++)
    {
        // step 1 take a half step in momentum
        ParameterDict grad = differ_(qs_);
        for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
            it->second -= epsilon_ / 2 * grad[it->first];

        // step 2 alternate full steps for position and momentum
        std::string reflectDim;
        double reflectTime;
        double t0 = 0;

        // reflectDim is set to "" if none exists
        FirstIntersectionTime(qs_, ps_, masses_, minima_, maxima_, epsilon_ - t0, reflectDim, reflectTime);
        while (reflectDim != "")
        {
            // use up the time to get to the boundary
            t0 += reflectTime;

            // get to this boundary
            for (ParameterDict::iterator it = qs_.begin(); it != qs_.end(); ++it)
                it->second += reflectTime * ps_[it->first] / masses_.at(it->first);

            // flip the offending momentum
            ps_.at(reflectDim) *= -1;

            // on to the next reflection if there is one
            FirstIntersectionTime(qs_, ps_, masses_, minima_, maxima_, epsilon_ - t0, reflectDim, reflectTime);
        }

        // finally advance the remaining time knowing there are no futher reflection

        double remainingTime = epsilon_ - t0;
        for (ParameterDict::iterator it = qs_.begin(); it != qs_.end(); ++it)
            it->second += remainingTime * ps_[it->first] / masses_.at(it->first);

        // step 3 make a half step for momentum
        grad = differ_(qs_);
        for (ParameterDict::iterator it = ps_.begin(); it != ps_.end(); ++it)
            it->second -= epsilon_ / 2 * grad[it->first];
    }
}
