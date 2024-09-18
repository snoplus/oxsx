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
class DenseMatrix;
class Convolution : public Systematic
{
public:
   Convolution(const std::string &name_) : fDist(NULL), fCachedPermutationMatrix(false), fName(name_), fPermMatrixIdentity(true) {}
   ~Convolution();
   void SetFunction(PDF *function_);
   void SetConditionalPDF(ConditionalPDF *function_);
   void Construct();

   // Make this fittable, by delegating to the underlying pdf
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
   bool fCachedPermutationMatrix;
   AxisCollection fSubMapAxes;
   AxisCollection fNotSubMapAxes;
   // std::vector<std::vector<size_t>> fCompatibleBins;
   // the systematic subMap bin for each global bin of pdf
   // std::vector<size_t> fSysBins;
   std::string fName;
   bool fPermMatrixIdentity;
   SparseMatrix fBinningPerm;
   SparseMatrix fBinningPermT;

   AxisCollection DetermineAxisSubCollection(const std::vector<size_t>& rel_indices);
   size_t BlockedBinningIndex(size_t bin_index, const std::vector<size_t>& relativeIndices);
   void CachePermutationMatrix();
   // void CacheCompatibleBins();
   void ConstructSubmatrix(std::vector<long long unsigned int>& column_indices, std::vector<long long unsigned int>& row_indices,
                           std::vector<double>& vals) const;
   void MakeBlockMatrix(const std::vector<long long unsigned int>& column_indices, const std::vector<long long unsigned int>& row_indices,
                           const std::vector<double>& vals, SparseMatrix& response_blocked);
};
#endif
