#ifndef __OXSX_TYPE_TRAITS__
#define __OXSX_TYPE_TRAITS__
#include <limits>

namespace TypeTraits
{

    // this is a back port of enable_if from c++11, implementation from c++ docs
    // substitution always works here, regardless of the value of B
    template <bool B, typename T = void>
    struct enable_if
    {
    };

    // this subsitution can override the previous one if B == true, then it will give it a typedef T
    template <typename T>
    struct enable_if<true, T>
    {
        typedef T type;
    };

    // same for conditional
    template <bool B, class T, class F>
    struct conditional
    {
        typedef T type;
    };

    template <class T, class F>
    struct conditional<false, T, F>
    {
        typedef F type;
    };

    // the next three are taken from https://gist.github.com/louisdx/1076849 - credit there
    // does this object have a const iterator
    template <typename T>
    struct has_const_iterator
    {
    private:
        typedef char one;
        typedef struct
        {
            char array[2];
        } two;

        template <typename C>
        static one test(typename C::const_iterator *);
        template <typename C>
        static two test(...);

    public:
        static const bool value = sizeof(test<T>(0)) == sizeof(one);
        typedef T type;
    };

    // does this object have a begin and an end
    template <typename T>
    struct has_begin_end
    {
        struct Dummy
        {
            typedef void const_iterator;
        };
        typedef typename conditional<has_const_iterator<T>::value, T, Dummy>::type TType;
        typedef typename TType::const_iterator iter;

        struct Fallback
        {
            iter begin() const;
            iter end() const;
        };
        struct Derived : TType, Fallback
        {
        };

        template <typename C, C>
        struct ChT;

        template <typename C>
        static char (&f(ChT<iter (Fallback::*)() const, &C::begin> *))[1];
        template <typename C>
        static char (&f(...))[2];
        template <typename C>
        static char (&g(ChT<iter (Fallback::*)() const, &C::end> *))[1];
        template <typename C>
        static char (&g(...))[2];

        static bool const beg_value = sizeof(f<Derived>(0)) == 2;
        static bool const end_value = sizeof(g<Derived>(0)) == 2;
    };

    template <typename T>
    struct is_container
    {
        static const bool value = has_const_iterator<T>::value &&
                                  has_begin_end<T>::beg_value && has_begin_end<T>::end_value;
    };

    // fundamental data type stuff, use the types for which numeric limits is specialized
    template <typename T>
    struct is_number
    {
        static const bool value = std::numeric_limits<T>::is_specialized;
    };

} // namespace
#endif
