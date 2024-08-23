// Used to parse .ini configuration files into c++ types
// Supports all basic types & strings and stl (non-associative) containers except forward_list
// if used on a container existing elements are not removed
// example use:

// filename: example.ini
//
// myint = 5
// mystr = hello

// filename: client.cpp
//
// #include <ini.hpp>
// #include <ConfigLoader.hh>
// #include <string>
// #include <iostream>
//
// ConfigLoader::Open("example.ini");
// ConfigLoader::Load("myint", myint);
// ConfigLoader::Load("myint", mystr);
// ConfigLoader::Close();
// std::cout << myint == 5 << std::endl; // 1
// std::cout << mystr == "hello" << std::endl; // 1

#ifndef __OXSX__CONFIG_LOADER__
#define __OXSX__CONFIG_LOADER__
#include <TypeTraits.hpp>
#include <ini.hpp>
#include <set>
#include <string>

using TypeTraits::enable_if;
using TypeTraits::is_container;
using TypeTraits::is_number;

class ConfigLoader
{
public:
  static void Open(const std::string &fileName_);
  static void Close();

  // Fundamental numeric types
  template <typename TargetType>
  static void
  Load(const std::string &section_, const std::string &fieldName_, TargetType &loadVal_, typename enable_if<is_number<TargetType>::value, int>::type = 0);

  // another one for containers
  template <typename TargetType>
  static void
  Load(const std::string &section_, const std::string &fieldName_, TargetType &loadVal_, typename enable_if<is_container<TargetType>::value, int>::type = 0);

  //
  static void
  Load(const std::string &section_, const std::string &fieldName_, std::string &);

  static std::set<std::string> ListSections();

private:
  static void CheckExists(const std::string &section_, const std::string &field_);
  static INI::Parser *fParser;
};

// nested class performs string conversions
template <typename TargetType>
struct Converter
{
  TargetType operator()(const std::string &s_) const;
};

// specialise this for strings
template <>
struct Converter<std::string>
{
  std::string operator()(const std::string &s_) const { return s_; }
};

// convert a whole container
template <typename InContainer, typename OutContainer, typename ConverterSp>
static void
ConvertContainer(const InContainer &incntr_, OutContainer &outcntr_, const ConverterSp &cnvtr_);

#include <ConfigLoader.hpp>

#endif
