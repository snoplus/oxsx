#ifndef __OXSX_EVENT_CONVOLUTION__
#define __OXSX_EVENT_CONVOLUTION__
#include <EventSystematic.h>
#include <vector>
#include <string>
#include <iostream>

class PDF;
class ConditionalPDF;
class Event;
class EventConvolution : public EventSystematic
{
public:
  // Constructors
  EventConvolution(const std::string &name_) : fDist(NULL), fName(name_) {}
  ~EventConvolution();
  EventConvolution(const EventConvolution &);
  EventConvolution operator=(const EventConvolution &);

  // Get/Set
  void SetPDF(PDF *f);
  void SetConditionalPDF(ConditionalPDF *f);

  // Event Systematic interface
  Event operator()(const Event &);

  // Fit Component interface defers to underlying CPDF
  void SetParameter(const std::string &name_, double value);
  double GetParameter(const std::string &name_) const;

  void SetParameters(const ParameterDict &);
  ParameterDict GetParameters() const;
  size_t GetParameterCount() const;

  std::set<std::string> GetParameterNames() const;
  void RenameParameter(const std::string &old_, const std::string &new_);

  std::string GetName() const;
  void SetName(const std::string &);

private:
  ConditionalPDF *fDist;
  std::string fName;
};
#endif
