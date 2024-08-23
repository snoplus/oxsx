#ifndef __OXSX_COMBINATIONS__
#define __OXSX_COMBINATIONS__
#include <Exceptions.h>
#include <vector>
#include <utility>
#include <assert.h>
#include <set>

namespace Combinations
{

    template <typename T1>
    std::set<std::pair<typename T1::value_type, typename T1::value_type>>
    AllCombsNoDiag(const T1 &v1_)
    {

        typedef typename T1::value_type ValueType;
        typedef typename std::set<std::pair<ValueType, ValueType>> SetType;
        typedef typename T1::const_iterator It;

        SetType pairs;
        for (It it1 = v1_.begin(); it1 != v1_.end(); ++it1)
        {
            for (It it2 = v1_.begin(); it2 != it1; ++it2)
            {
                pairs.insert(make_pair(*it1, *it2));
            }
        }
        return pairs;
    }

    template <typename T1>
    std::vector<T1> Range(const T1 &N, const T1 &start = 0)
    {
        assert(N > start);

        std::vector<T1> ret;
        ret.reserve(N - start);
        for (T1 i = start; i < N; i++)
            ret.push_back(i);
        return ret;
    }

    void RecursiveCombinations(size_t depth, std::vector<std::vector<size_t>> &combs, std::vector<size_t> vec, std::vector<size_t> maxes);

}
#endif
