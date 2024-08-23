#include <Combinations.hpp>
#include <Exceptions.h>
#include <vector>
#include <utility>
#include <assert.h>
#include <set>
#include <iostream>

void Combinations::RecursiveCombinations(size_t depth, std::vector<std::vector<size_t>> &combs, std::vector<size_t> vec, std::vector<size_t> maxes)
{

  if (depth > 0)
  {
    for (size_t i = 0; i < maxes[depth - 1]; i++)
    {
      vec[depth - 1] = i;
      RecursiveCombinations(depth - 1, combs, vec, maxes);
    }
  }
  else
  {
    combs.push_back(vec);
  }
  return;
}
