/*****************************************/
/* A simple shift error on an observable */
/*****************************************/
#ifndef __OXSX_SHIFT__
#define __OXSX_SHIFT__
#include <Systematic.h>
#include <string>
#include <DenseMatrix.h>

class Shift : public Systematic
{
public:
   Shift(const std::string &name_) : fShift(1), fName(name_), fParamName("shift"), fCachedBinMapping(false) {}
   void SetShift(double);
   double GetShift() const;

   void Construct();

   // Adjustable shift
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
   double fShift;
   std::string fName;
   std::string fParamName;

   std::vector<size_t> fDistTransBinMapping;
   DenseMatrix fMappingDistAndTrans;
   bool fCachedBinMapping;

   void CacheBinMapping();
   double GetBinContribution(const BinAxis &shiftAxis, double shiftedLow,
                             double shiftedHigh, size_t bin_test);
};
#endif
