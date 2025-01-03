/*****************************************************************************************************/
/* Convolution of the pdf with another compatible pdf. The pdf must have dimensionality smaller than */
/* the pdfs it will act on, and must implement IntegrablePdf method Integrate() to allow the         */
/* convolution to integrate it over bin boundaries.                                                  */
/* The parameters set in this systematic are forwarded on to the internal pdf, you then need to call */
/* Construct() to update the response matrix                                                         */
/*****************************************************************************************************/

#ifndef __OXSX_CONVOLUTION__
#define __OXSX_CONVOLUTION__
#include <Systematic.h>

class ConditionalPDF;
class PDF;

class Convolution : public Systematic
{
public:
   // Constructor/Destructor
   Convolution(const std::string &name_) : fDist(NULL), fName(name_), fCachedIndexPermutations(false) {}
   ~Convolution();
   // Two ways of setting the kernel for convolution:
   void SetFunction(PDF *function_);                  // Just give me a PDF object directly...
   void SetConditionalPDF(ConditionalPDF *function_); // ...Or create a ConditionalPDF object first yourself
   // Construct the response matrix for this convolution object! (Systematic Interface)
   void Construct();

   // FitComponent interface
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
   ConditionalPDF *fDist; // kernel used in convolution
   std::string fName;     // name of this object
   // Members needed for efficient calculations:
   bool fCachedIndexPermutations;                    // have we run the CacheIndexPermutations() method yet?
   AxisCollection fSubMapAxes;                       // axes in fAxes associated with fTransObs
   AxisCollection fNotSubMapAxes;                    // axes in fAxes NOT associated with fTransObs
   std::vector<long long unsigned int> fColumnPerms; // Mapping from blocked matrix bin index -> true bin index

   AxisCollection DetermineAxisSubCollection(const std::vector<size_t> &rel_indices) const;
   size_t BlockedBinningIndex(size_t bin_index, const std::vector<size_t> &relativeIndices) const;
   void CacheIndexPermutations();
   void ConstructSubmatrix(std::vector<long long unsigned int> &column_indices, std::vector<long long unsigned int> &row_indices,
                           std::vector<double> &vals) const;
   void MakeFullMatrix(const std::vector<long long unsigned int> &column_indices, const std::vector<long long unsigned int> &row_indices,
                       const std::vector<double> &vals, SparseMatrix &response_blocked) const;
};
#endif
