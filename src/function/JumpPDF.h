#ifndef __OXSX_JUMP_PDF__
#define __OXSX_JUMP_PDF__
#include <ConditionalPDF.h>
#include <vector>
#include <string>

class PDF;
class JumpPDF : public ConditionalPDF{
 public:  
  // Constructors/Destructor
  JumpPDF(const std::string& name_): fPDF(NULL), fName(name_) {}
  JumpPDF(const std::string& name_, PDF*);
  ~JumpPDF();
  JumpPDF(const JumpPDF&);
  JumpPDF operator=(const JumpPDF&);
  ConditionalPDF* Clone() const;

  // ConditDist Interface, just execute the function on x - x2
  virtual double ConditionalProbability(const std::vector<double>& x_,
                                        const std::vector<double>& x2_);
  
  virtual double Integral(const std::vector<double>& mins_,
                          const std::vector<double>& maxs_,
                          const std::vector<double>& x2_) const;
  
  std::vector<double> Sample(const std::vector<double>& x2_) const;

  // Elementwise difference/sum between x and x2
  std::vector<double> Diff(const std::vector<double>& x_, 
                           const std::vector<double>& x2_) const;

  std::vector<double> Sum(const std::vector<double>& x_, 
                          const std::vector<double>& x2_) const;

  // Fit Component Interface - defer to the underlying unconditional PDF
  void   SetParameter(const std::string& name_, double value);
  double GetParameter(const std::string& name_) const;
  
  void   SetParameters(const ParameterDict&);
  ParameterDict GetParameters() const;
  size_t GetParameterCount() const;
  
  std::set<std::string> GetParameterNames() const;
  void   RenameParameter(const std::string& old_, const std::string& new_);
  
  std::string GetName() const;
  void SetName(const std::string&);
  
 private:
  PDF* fPDF;
  std::string fName;
};
#endif
