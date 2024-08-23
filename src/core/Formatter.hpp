/*
A wrapper around stringstream. stringstream << returns a ostream pointer,
that doesn't allow implicit conversion to strings for error
message formatting - this does. See stack overflow 12261915
*/
#ifndef __OXSX_FORMATTER__
#define __OXSX_FORMATTER__
#include <sstream>

class Formatter
{
public:
    Formatter() {}
    ~Formatter() {}
    template <typename Type>
    Formatter &operator<<(const Type &val_)
    {
        fStringStream << val_;
        return *this;
    }
    std::string str() const { return fStringStream.str(); }
    operator std::string() const { return fStringStream.str(); }

private:
    std::stringstream fStringStream;
    Formatter(const Formatter &);
    Formatter operator=(const Formatter &);
};
#endif
