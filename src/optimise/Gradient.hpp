#include <iostream>

template <typename Stat>
ParameterDict
Gradient<Stat>::operator()(const ParameterDict &x0_)
{
    ParameterDict grad;

    for (ParameterDict::const_iterator it = x0_.begin(); it != x0_.end(); ++it)
    {
        ParameterDict high = x0_;
        high[it->first] += fEpsilon;

        fStat.SetParameters(high);
        double highVal = fStat.Evaluate();

        ParameterDict low = x0_;
        low[it->first] -= fEpsilon;

        fStat.SetParameters(low);
        double lowVal = fStat.Evaluate();

        grad[it->first] = (highVal - lowVal) / 2 / fEpsilon;
    }

    return grad;
}
