/**************************************************************************************************/
/* A BinnedPdf is an ND histogram along with a data representation that tells the pdf where it    */
/* looks for the data to pass to the histogram inside. Most methods just defer to the histogram.  */
/* The two are serparate classes, because other things are histograms e.g. LH space but not pdfs  */
/* in the sense that they have a data representation and examine events                           */
/**************************************************************************************************/

#ifndef __OXSX_BINNED_ED__
#define __OXSX_BINNED_ED__
#include <Histogram.h>
#include <EventDistribution.h>
#include <vector>
#include <map>
#include <string>
#include <limits>

class Event;
class BinnedED : public EventDistribution
{
public:
   BinnedED() {}
   BinnedED(const std::string &name_, const AxisCollection &axes_);
   BinnedED(const std::string &name_, const Histogram &histo_);
   EventDistribution *Clone() const;

   double Probability(const Event &) const;
   double Probability(const std::vector<double> &) const;
   double Integral() const;
   void Normalise();
   void Scale(double s_);

   void Fill(const std::vector<double> &vals_, double weight_ = 1);
   void Fill(const Event &data_, double weight_ = 1);
   void Fill(double val_, double weight_ = 1);

   size_t FindBin(const std::vector<double> &vals_) const;
   size_t FindBin(const Event &data_) const;

   std::vector<size_t> UnpackIndices(size_t bin_) const;
   size_t FlattenIndices(const std::vector<size_t> &indices_) const;

   const AxisCollection &GetAxes() const;
   void SetAxes(const AxisCollection &axes_);

   std::vector<double> GetBinContents() const;
   void SetBinContents(const std::vector<double> &data_);

   std::vector<double> Means() const;
   std::vector<double> Variances() const;

   const Histogram &GetHistogram() const;
   void SetHistogram(const Histogram &);

   double GetBinContent(size_t bin_) const;
   size_t GetNBins() const;

   void AddBinContent(size_t bin_, double content_);
   void SetBinContent(size_t bin_, double content_);
   void Empty();
   unsigned GetNDims() const;

   BinnedED Marginalise(const std::vector<std::string> &indices_) const;
   BinnedED Marginalise(const std::string &index_) const;

   BinnedED GetSlice(const std::map<std::string, size_t> &namesAndBins_) const;

   void SetObservables(const std::vector<std::string> &);
   const std::vector<std::string> &GetObservables() const;

   void Add(const BinnedED &other_, double weight_ = 1);
   void Multiply(const BinnedED &other_);
   void Divide(const BinnedED &other_);

   std::string GetName() const;
   void SetName(const std::string &);

   void AddPadding(double padding_ = std::numeric_limits<double>::min());

protected:
   ObsSet fObservables;
   Histogram fHistogram;
   std::string fName;
};
#endif
