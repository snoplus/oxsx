/*****************************************/
/* A scale function error on an observable */
/*****************************************/
#ifndef __OXSX_SCALEFUNCTION__
#define __OXSX_SCALEFUNCTION__
#include <Systematic.h>
#include <string>
#include <DenseMatrix.h>
#include <functional>

/*
 * Firstly, create a typdef for a general function that just
 * takes in a double and the parameter dict, and returns a
 * double.
 */
typedef std::function<double(const ParameterDict &,
                             const double &)>
    ScaleFunc;

class ScaleFunction : public Systematic
{
public:
   ScaleFunction(const std::string &name_) : fScaleFunc(), fParamDict(), fName(name_), fCachedBinMapping(false) {}

   void SetScaleFunction(const ScaleFunc &scale_func,
                         const std::vector<std::string> &param_names);

   void Construct();

   void SetParameter(const std::string &name_, double value);
   double GetParameter(const std::string &name_) const;

   void SetParameters(const ParameterDict &);
   ParameterDict GetParameters() const { return fParamDict; }
   size_t GetParameterCount() const { return fParamDict.size(); }

   std::set<std::string> GetParameterNames() const;
   void RenameParameter(const std::string &old_, const std::string &new_);

   std::string GetName() const { return fName; }
   void SetName(const std::string &name_) { fName = name_; }

private:
   ScaleFunc fScaleFunc;
   ParameterDict fParamDict;
   std::string fName;
   std::string fParamName;

   std::vector<size_t> fDistTransBinMapping;
   DenseMatrix fMappingDistAndTrans;
   bool fCachedBinMapping;

   void CacheBinMapping();
   double GetBinContribution(const BinAxis &scaleAxis, double scaledLow,
                             double scaledHigh, size_t bin_test);
};
#endif
