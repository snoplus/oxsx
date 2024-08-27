#ifndef __OXSX_COMPARISION__
#define __OXSX_COMPARISION__
template <typename T1, typename T2, typename Comp = std::less<T2>>
struct PairSortBySecond
{
  bool operator()(const std::pair<T1, T2> &left_,
                  const std::pair<T1, T2> &right_) const
  {
    Comp c;
    return c(left_.second, right_.second);
  }
};
#endif
