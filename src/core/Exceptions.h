#ifndef __OXSX_EXCEPTIONS__
#define __OXSX_EXCEPTIONS__
#include <stdexcept>
#include <Formatter.hpp>

class OXSXException : public std::runtime_error
{
public:
    OXSXException(const std::string &errorStr_) : runtime_error(errorStr_) {}
};

// Wrong dimensionality
class DimensionError : public OXSXException
{
public:
    DimensionError(const std::string &errorStr_) : OXSXException(errorStr_) {}
    DimensionError(const std::string &caller_, int expected_, int recieved_,
                   const std::string &errorStr_ = "")
        : OXSXException(
              Formatter() << caller_
                          << " Expected: " << expected_
                          << " Got: " << recieved_
                          << " " << errorStr_) {}
};

// Invalid value
class ValueError : public OXSXException
{
public:
    ValueError(const std::string &errorStr_) : OXSXException(errorStr_) {}
};

// missing from container
class NotFoundError : public OXSXException
{
public:
    NotFoundError(const std::string &errorStr_) : OXSXException(errorStr_) {}
};

// Failure in i/o operations
class IOError : public OXSXException
{
public:
    IOError(const std::string &errorStr_) : OXSXException(errorStr_) {}
};

// failure in ObsSet system
class RepresentationError : public OXSXException
{
public:
    RepresentationError(const std::string &errorStr_) : OXSXException(errorStr_) {}
};

// failure in fit parameter structure
class ParameterError : public OXSXException
{
public:
    ParameterError(const std::string &errorStr_) : OXSXException(errorStr_) {}
    ParameterError(const std::string &className_,
                   int expected_,
                   int recieved_,
                   const std::string &extraMessage_ = "")
        : OXSXException(Formatter() << className_ << " Expected: " << expected_ << " Got: " << recieved_ << " "
                                    << extraMessage_) {}
};

// Initialisation Error
class LogicError : public OXSXException
{
public:
    LogicError(const std::string &errorStr_) : OXSXException(errorStr_) {}
};

// NULL pointer access
class NULLPointerAccessError : public OXSXException
{
public:
    NULLPointerAccessError(const std::string &errorStr_) : OXSXException(errorStr_) {}
    NULLPointerAccessError(const std::string &caller_,
                           const std::string &errorStr_)
        : OXSXException(
              Formatter() << caller_
                          << " Attempted access on NULL pointer. " << errorStr_) {}
};

class ConfigSectionMissing : public OXSXException
{
public:
    ConfigSectionMissing(const std::string &str_) : OXSXException(str_) {}
};

class ConfigFieldMissing : public OXSXException
{
public:
    ConfigFieldMissing(const std::string &str_) : OXSXException(str_) {}
};

#endif
