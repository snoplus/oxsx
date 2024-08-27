#include <stdexcept>
//////////////////////////
// CLASS IMPLEMENTATION //
//////////////////////////

template <typename Container>
double
ContainerParameter<Container>::Get() const
{
    if (fIndex >= fContainer.size())
        throw std::out_of_range("Container Parameter is out of bounds!");
    return GetNth(fContainer, fIndex, Category());
}

template <typename Container>
void ContainerParameter<Container>::Set(double val_)
{
    if (fIndex >= fContainer.size())
        throw std::out_of_range("Container Parameter is out of bounds!");
    GetNth(fContainer, fIndex, Category()) = val_;
}

//////////////////////////////////////////////////////////////////
// STANDALONE FUNCTION IMPLEMENTIONS OF GetNTth, compiler will  //
// choose one for ContainerParameter based on the category tag  //
// for the container iterator                                   //
//////////////////////////////////////////////////////////////////

// if you have random access just use it!
template <typename Container>
double &
GetNth(Container &C_, size_t n_, std::random_access_iterator_tag)
{
    return C_[n_];
}

// otherwise we have only forward and backward motion
template <typename Container>
double &
GetNth(Container &C_, size_t n_, std::bidirectional_iterator_tag)
{
    typename Container::iterator it = C_.begin();
    for (size_t i = 0; i < n_; i++)
        ++it;
    return *it;
}

// or just forward motion
template <typename Container>
double &
GetNth(Container &C_, size_t n_, std::forward_iterator_tag)
{
    return GetNth(C_, n_, std::bidirectional_iterator_tag());
}
