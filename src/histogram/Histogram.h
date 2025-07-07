/*************************************************************************************/
/* An ND histogram is a combination of an axis collection that defines the binning   */
/* and a single vector of bin contents, indexed by the corresponding global bin ID.  */
/* Conversions between global bin ID and the indices on each axis are performed by  */
/* the axis collection                                                               */
/*************************************************************************************/
#ifndef __OXSX_HISTOGRAM__
#define __OXSX_HISTOGRAM__
#include <AxisCollection.h>
#include <vector>
#include <map>
#include <limits>

class Histogram
{
public:
   Histogram() : fNBins(0), fNDims(0) {}
   Histogram(const AxisCollection &axes_);

   double Integral() const;
   double Integral(double,double);//calculates the Bounded Integral for 1-D Histograms only
   void Normalise();
   void Scale(double);

   void Fill(const std::vector<double> &vals_, double weight_ = 1);
   void Fill(const std::map<std::string, double> &vals_, double weight_ = 1);
   void Fill(double val_, double weight_ = 1);

   size_t FindBin(const std::vector<double> &vals_) const;
   size_t FindBin(double);
  
   std::vector<size_t> UnpackIndices(size_t bin_) const;
   size_t FlattenIndices(const std::vector<size_t> &indices_) const;

   const AxisCollection &GetAxes() const;
   void SetAxes(const AxisCollection &axes_);

   std::vector<double> GetBinContents() const;
   void SetBinContents(const std::vector<double> &data_);

   std::vector<double> Means() const;
   std::vector<double> Variances() const;

   double GetBinContent(size_t bin_) const;
   size_t GetNBins() const;
   size_t GetNDims() const;
   void AddBinContent(size_t bin_, double content_);
   void SetBinContent(size_t bin_, double content_);
   void Empty();

   Histogram Marginalise(const std::vector<std::string> &indices_) const;
   Histogram Marginalise(const std::string &index_) const;

   Histogram GetSlice(const std::map<std::string, size_t> &fixedBins_) const;

   double GetBinLowEdge(size_t bin_, size_t dim_) const;
   double GetBinHighEdge(size_t bin_, size_t dim_) const;
   double GetBinCentre(size_t bin_, size_t dim_) const;

   void Add(const Histogram &, double weight = 1);
   void Multiply(const Histogram &);
   void Divide(const Histogram &);

   std::vector<std::string> GetAxisNames() const;
   size_t GetAxisIndex(const std::string &name_) const;
   void AddPadding(double padding_ = std::numeric_limits<double>::min());

private:
   AxisCollection fAxes;
   std::vector<double> fBinContents;
   size_t fNBins;
   size_t fNDims;
   void Recurse(size_t numFreeIdx, std::vector<size_t> binsEachAxis, std::vector<size_t> coords, std::vector<std::vector<size_t>> &localIdx, std::vector<size_t> free_ax) const;
};
#endif
